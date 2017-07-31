""" This module contains functions to be use for Phonon calculations.
  
  get_spacegroup:
    define an equivalent to spglib.get_spacegroup
    Try out with:
        from ase.lattice import bulk
        import get_spacegroup
        atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
        sg = ifit.get_spacegroup(atoms)
        sg
            Spacegroup(225, setting=1)
        atoms.info
            {'spacegroup': Spacegroup(225, setting=1)}
            
  phonons_run:
    same as Phonons.run but using the spacegroup information to speed-up the
    calculation.
"""

# we need spacegroup from ASE
try:
    from ase.spacegroup import Spacegroup         # For ASE version 3.10 or later
except ImportError:
    from ase.lattice.spacegroup import Spacegroup # For ASE version 3.9 or before

import sys
import pickle
import numpy
import warnings
import os

from ase.atoms     import Atoms
from ase.utils     import opencew
from ase.parallel  import rank, barrier

# check if we have access to get_spacegroup from spglib
# https://atztogo.github.io/spglib/
has_spglib = False
try:
    import spglib                   # For version 1.9 or later
    has_spglib = True
except ImportError:
    try:
        from pyspglib import spglib # For versions 1.8.x or before
        has_spglib = True
    except ImportError:
        pass

# test if PhonoPy is there (which usually includes spglib)
has_phonopy = False
try:
    from phonopy import Phonopy
    has_phonopy = True
except ImportError:
    pass

# ------------------------------------------------------------------------------
def get_spacegroup(atoms, symprec=1e-5, method='phonopy'):
    """Determine the spacegroup to which belongs the Atoms object.
    
    Parameters:
    atoms:    an Atoms object
    symprec:  Symmetry tolerance, i.e. distance tolerance in Cartesian 
              coordinates to find crystal symmetry.
    method:   'phonopy' when available, or 'ase'
              
    The Spacegroup object is returned, and stored in atoms.info['spacegroup'] 
    when this key does not exist (avoids overwrite). To force overwrite of the 
    spacegroup first use:
        del atoms.info["spacegroup"]
                  
    Examples:
    
    >>> from ase.lattice import bulk
    >>> atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    >>> sg = ifit.get_spacegroup(atoms)
    """

    # use spglib when it is available (and return)
    if has_spglib and method is 'phonopy':
        sg    = spglib.get_spacegroup(atoms)
        sg_no = int(sg[sg.find("(")+1:sg.find(")")])
        atoms.info["spacegroup"] = Spacegroup(sg_no)
        return atoms.info["spacegroup"]
    
    # no spglib, we use our own spacegroup finder. Not as robust as spglib.
    # we center the Atoms positions on each atom in the cell, and find the 
    # spacegroup of highest symmetry
    found = None
    for kind, pos in enumerate(atoms.get_scaled_positions()):
        sg = _get_spacegroup(atoms, symprec=1e-5, center=kind)
        if found is None or sg.no > found.no:
            found = sg
    
    # return None when no space group is found (would be surprising)
    if found is not None and 'spacegroup' not in atoms.info:
          atoms.info["spacegroup"] = found

    return found
            
# ------------------------------------------------------------------------------
def _get_spacegroup(atoms, symprec=1e-5, center=None):
    """ASE implementation of get_spacegroup.
    """ 
    
    # we try all available spacegroups from 1 to 230, backwards
    # a Space group is the collection of all symmetry operations which lets the 
    # unit cell invariant.
    found      = None
    positions  = atoms.get_scaled_positions(wrap=True)  # in the lattice frame
    
    # make sure we are insensitive to translation. this choice is arbitrary and 
    # could lead to a 'slightly' wrong guess for the Space group, e.g. do not  
    # guess centro-symmetry.
    if center:
        try:
            positions -= positions[center]
        except IndexError:
            pass
    
    # search space groups from the highest symmetry to the lowest
    # retain the first match
    for nb in range(230,0,-1):
        sg        = Spacegroup(nb)
        #
        # now we scan all atoms in the cell and look for equivalent sites
        sites,kinds = sg.equivalent_sites(positions, 
                onduplicates='keep', symprec=symprec)
        #    
        # the equivalent sites should match all other atom locations in the cell
        # as the spacegroup transforms the unit cell in itself
        if len(sites) == len(positions):
            # store the space group into the list
            found = sg
            break

    return found
        
# ------------------------------------------------------------------------------
def _get_wigner_seitz(atoms, move=0, wrap=1, symprec=1e-6):
    """Compute the Wigner-Seitz cell.
       this routine is rather slow, but it allows to limit the number of moves
       
       code outrageously copied from PHON/src/get_wig by D. Alfe
       
       This is equivalent to a Voronoi cell, but centered on atom 'move'
       An alternative is to use scipy.spatial.Voronoi, but it does not properly
       center the cell. Slow but secure...
       
       input:
          atoms: an ASE Atoms object, supercell
          move:  the index  of the atom in the cell to be used as center
          wrap:  the extent of the supercell
       output:
          ws:    atom positions in the Wigner-Seitz cell, Cartesian coordinates
    """
    
    # this function is fully equivalent to the PHON/src/get_wig.f
    
    # get the cell definition
    L     = atoms.get_cell()
    # get fractional coordinates in the cell/supercell
    xtmp  = atoms.get_scaled_positions().copy()
    # set the origin on the 'moving' atom
    xtmp -= xtmp[move]

    # construct the WS cell (defined as the closest vectors to the origin)
    r = range(-wrap,wrap+1)
    for na in range(len(xtmp)):
        # convert from fractional to Cartesian
        xtmp[na] = numpy.dot(L.T, xtmp[na].T)
        temp1 = xtmp[na]
        for i in r:
            for j in r:
                for k in r:
                    temp = xtmp[na] + i*L[0] + j*L[1] + k*L[2]
                    if numpy.all(abs(numpy.linalg.norm(temp) < numpy.linalg.norm(temp1))):
                        temp1 = temp
        xtmp[na] = temp1
    
    return xtmp
    
# ------------------------------------------------------------------------------
def phonons_run_eq(phonon, supercell):
    """Run the calculation for the equilibrium lattice
    
    input:
    
      phonon: ASE Phonon object with Atoms and Calculator
      supercell: an ASE Atoms object, supercell
    
    return:
      output: forces (ndarray)
    """
    # Do calculation on equilibrium structure
    filename = phonon.name + '.eq.pckl'

    fd = opencew(filename)
    if fd is not None:
        # Call derived class implementation of __call__
        output = phonon.__call__(supercell)
                    
        # Write output to file
        if rank == 0:
            pickle.dump(output, fd, protocol=2)
            sys.stdout.write('Writing %s\n' % filename)
            fd.close()
            # check forces
            try:
                fmax = output.max()
                fmin = output.min()
                sys.stdout.write('[ASE] Equilibrium forces min=%g max=%g\n' % (fmin, fmax))
            except AttributeError:
                sys.stdout.write('[ASE] output is has no min/max (list)\n');
                pass
        sys.stdout.flush()
    else:
        # read previous data
        output = pickle.load(open(filename))
        
    return output
        
# ------------------------------------------------------------------------------
def phonons_run(phonon, single=True, difference='central'):
    """Run the calculations for the required displacements.
    This function uses the Spacegroup in order to speed-up the calculation.
    Every iteration is checked against symmetry operators. In case other moves 
    are imaged from a symmetry, the pickle files are written in advance, which 
    then allows to skip subsequent steps.

    This will do a calculation for 6 displacements per atom, +-x, +-y, and
    +-z. Only those calculations that are not already done will be
    started. Be aware that an interrupted calculation may produce an empty
    file (ending with .pckl), which must be deleted before restarting the
    job. Otherwise the calculation for that displacement will not be done.
    
    This implementation is the same as ASE, but allows to select the type of
    force gradient to use (difference). Also, it does make use of the spacegroup 
    to lower the number of displacements, using the strategy adopted in PHON 
    (D. Alfe). It is not as good as PhonoPy, but remains much simpler.
    
    input:
    
      phonon: ASE Phonon object with Atoms and Calculator
      single: when True, the forces are computed only for a single step, and then
              exit. This allows to split the loop in independent iterations. When
              calling again the 'run' method, already computed steps are ignored,
              missing steps are completed, until no more are needed. When set to
              False, all steps are done in a row.
      difference: the method to use for the force difference (gradient) computation.
              can be 'central','forward','backward'. The central difference is 
              more accurate but requires twice as many force calculations.
    
    output:
    
      True when a calculation step was performed, False otherwise or no more is needed.

    """
    
    # prepare to use symmetry operations
    # check if a spacegroup is defined, else find it
    if 'spacegroup' not in phonon.atoms.info or phonon.atoms.info["spacegroup"] is None:
        sg    = get_spacegroup(phonon.atoms)

    if difference == 'backward':
        signs = [-1]    # only compute one side, backward difference
    elif difference == 'forward':
        signs = [1]     # only compute one side, forward difference
    else:
        signs = [-1,1]   # use central difference

    # Atoms in the supercell -- repeated in the lattice vector directions
    # beginning with the last
    supercell = phonon.atoms * phonon.N_c
    
    # Set calculator if provided
    assert phonon.calc is not None, "Provide calculator in Phonon __init__ method"
    supercell.set_calculator(phonon.calc)
    
    # when not central difference, we check if the equilibrium forces are small
    # and will use the '0' forces in gradient
    if len(signs) == 1 and not os.path.isfile(phonon.name + '.eq.pckl'): 
        # Do calculation on equilibrium structure
        if rank == 0:
            print "[ASE] Computing equilibrium"
        phonons_run_eq(phonon, supercell)
        if single:
            return True # and some more iterations may be required

    # Positions of atoms to be displaced in the reference cell
    natoms = len(phonon.atoms)
    offset = natoms * phonon.offset
    pos    = supercell.positions[offset: offset + natoms].copy()
    pos0   = supercell.positions.copy()
    
    L    = supercell.get_cell()         # latticecell       = at
    invL = numpy.linalg.inv(L)      # no 2*pi multiplier here = inv(at) = bg.T
    
    for sign in signs:  # +-
        # guess displacements and rotations
        dxlist, rotlist = _phonons_get_displacements(phonon, sign)
        
        for a in phonon.indices:    # atom index to move in cell
        
            # we determine the Wigner-Seitz cell around atom 'a' (Cartesian)
            ws = _get_wigner_seitz(supercell, a, 1)
        
            # check if this step has been done before (3 files written per step)
            filename = '%s.%d%s%s.pckl' % \
                       (phonon.name, a, 'xyz'[0], ' +-'[sign])
            if os.path.isfile(filename):
                    # Skip if these 3 moves are already done.
                    continue
            
            # compute the forces for 3 independent moves
            force0 = None
            force1 = []
            for index in range(len(dxlist)):
                # we first determine a set of independent displacements, using 
                # as much as possible the symmetry operators.
                disp = dxlist[index]
                rot  = rotlist[index]
                
                if force0 is not None and rot is not None:  # re-use previous forces
                    # we will now apply the rotation to the force array
                    output = _phonons_run_force1_rot(force0, ws, invL, rot, symprec=1e-6)
                    if output is None:
                        # failed using symmetry to derive force. Trigger full computation.
                        force0 = None
                    elif rank == 0:
                        print "[ASE] Imaging atom #%-3i %-3s    to " % \
                            (offset + a, supercell.get_chemical_symbols()[a]), pos[a] + disp, \
                            " (Angs) using rotation:"
                        print rot
                        
                if force0 is None or rot is None: # compute forces
                    # move atom 'a' by 'disp'
                    supercell.positions[offset + a] = pos[a] + disp
                    if rank == 0:
                        print "[ASE] Moving  atom #%-3i %-3s    to " % \
                            (offset + a, supercell.get_chemical_symbols()[a]), pos[a] + disp, " (Angs)"
                        
                    # Call derived class implementation of __call__
                    output = phonon.__call__(supercell)
                    
                    # Return to initial positions
                    supercell.positions[offset + a] = pos[a]
                    
                    # store the forces for subsequent moves using symmetry
                    force0 = output
                
                # append the forces to the force1 list
                if output is None:
                    print "[ASE] Warning: force1 is None !!"

                force1.append(output)
                
            # when exiting the for 'i' loop, we have 3 independent 'force1' array
            # derive a Cartesian basis, and write pickle files
            force2 = _phonons_run_force2(phonon, dxlist, force1, a, sign, symprec=1e-6)
            
            # then we derive the Cartesian basis 'force2' array and write files
            if single:
                return True # and some more iterations may be required

    return False  # nothing left to do

# ------------------------------------------------------------------------------
def _phonons_move_is_independent(dxlist, dx, symprec=1e-6):
    """Test if a vector is independent compared to those in a list.
       The test is for collinearity, and singularity of the formed basis
       
       input:
           dxlist:  a list of vectors
           dx:      a new vector which is tested
           symprec: precision for comparing vectors
       output:
           True when dx is independent from dxlist, False otherwise.
    """

    # test if the rotated displacement is collinear to
    # a stored one (in list 'dxlist'). test is done on normalised vectors.
    if dx is None:
        return False
    
    dxnorm      = numpy.asarray([x/numpy.linalg.norm(x) for x in dxlist])
    iscollinear = False
    for index in range(len(dxnorm)):
        if numpy.linalg.norm(numpy.cross(dx/numpy.linalg.norm(dx), dxnorm[index])) < symprec:
            iscollinear = True
            break
    if iscollinear:
        return False  # collinear
            
    # Test for singular matrix when adding new move
    dxlist2 = dxlist[:] # copy the list
    dxlist2.append(dx)
    dxnorm = numpy.asarray([x/numpy.linalg.norm(x) for x in dxlist2])
    if len(dxlist2) == 3 and numpy.abs(numpy.linalg.det(dxnorm)) < symprec:
        return False # singular
        
    return True
    
# ------------------------------------------------------------------------------
def _phonons_rotated_displacements(disp, sg):
    """Check if a displacement can be imaged into other ones using the symmetry
       operators of the spacegroup.
       
       input:
           disp: an initial displacement (3-vector)
           sg:   ASE Space group
       output:
           dxlist:  list of independent displacements imaged from 'disp'
           rotlist: list of corresponding rotations from 'disp'
    """

    # we try all symmetry operations in the spacegroup
    dxlist = [disp]
    rotlist= [None]
    for rot, trans in sg.get_symop():
    
        # find the equivalent displacement from initial displacement
        # and rotation. First 'rot' is identity.
        dx = numpy.dot(rot, disp) # in Cartesian: dx = rot * disp
        
        # check if that rotated move contributes to an independent basis set
        if _phonons_move_is_independent(dxlist, dx):
            # store dx, rot; not work for disp itself that was added before
            dxlist.append(dx)
            rotlist.append(rot)
        
        # exit when 3 moves have been found from same generator move
        if len(dxlist) == 3:
            break
            
    return dxlist, rotlist   

# ------------------------------------------------------------------------------
def _phonons_get_displacements(phonon, sign):
    """Determine a set of displacements that best make use of crystal symmetries.
    
       input:
           phonon: ASE phonon containing an Atoms and a Spacegroup.
           sign:   -1 or +1, to indicate in which direction we move the atoms
       output:
           dxlist:  list of independent displacements
           rotlist: list of corresponding rotations from 'disp'
    """
    
    dxlist = [] # hold list of lists for tentative displacements and rotations
    rotlist= []
    L      = phonon.atoms.get_cell()
    sg     = phonon.atoms.info["spacegroup"]    # spacegroup for e.g. spglib
    
    # Determine displacement vectors to use (list[3]). 
    # We look for the equivalent displacements. We only add when symmetry 
    # operators provide equivalent sites (saves time).
    
    # we try with axes [xyz] or lattice axes that have most symmetries
    for i in range(6):
        if i < 3: # [xyz] directions
            disp     = numpy.asarray([0.0,0,0])
            disp[i] += sign * phonon.delta # in Angs, Cartesian
        else:     # lattice cell axes
            disp = L[i-3]/numpy.linalg.norm(L[i-3]) * sign * phonon.delta
        # test if symmetries can generate other independent moves
        this_dx, this_rot = _phonons_rotated_displacements(disp, sg)
        dxlist.append(this_dx)  # append lists
        rotlist.append(this_rot)
    
    # now we sort the lists by length of the sublists, i.e. select the moves 
    # that generate most equivalent moves by symmetry
             
    # get the index to use for sorting (decreasing order)
    lenlist = [len(x) for x in dxlist]
    order   = sorted(range(len(lenlist)), key=lambda k: lenlist[k], reverse=True)
    # reorder lists by decreasing size
    dxlist  = [ dxlist[j]  for j in order]
    rotlist = [ rotlist[j] for j in order]
    # now catenate all lists
    dxlist  = [j for i in dxlist  for j in i]
    rotlist = [j for i in rotlist for j in i]
    
    # and finally extract iteratively 3 vectors which are independent
    dxlist2 = []
    rotlist2= []
    for index in range(len(dxlist)):
        if _phonons_move_is_independent(dxlist2, dxlist[index]):
            dxlist2.append(dxlist[index])
            rotlist2.append(rotlist[index])
        if len(dxlist2) == 3:
            break
    
    # return only the first 3 independent entries
    return dxlist2, rotlist2

# ------------------------------------------------------------------------------
def _phonons_run_force1_rot(force, ws, invL, rot, symprec=1e-6):
    """From a given force set, we derive a rotated force set 
       by applying a single corresponding symmetry operator.
       
       code outrageously adapted from PHON/src/set_forces by D. Alfe
       
       input:
          force:    an array of Forces, size [natoms, 3]
          ws:       the Wigner-Seitz positions, supercell (Cartesian), obtained 
                    from _get_wigner_seitz(supercell)
          invL:     the inverse of the lattice supercell
          rot:      the rotation matrix
          symprec:  the precision for comparing positions
   
       output:
          nforce:   the rotated forces set
    """
    
    # this routine has been validated and behaves as in PHON/set_forces.f
    # the force is properly rotated when given the 'rot' matrix.
    
    # exit if no WS cell defined
    if ws is None:
        return None
        
    # we will now apply the rotation to the force array
    nforce = force.copy()

    # scan all atoms to compute the rotated force matrix
    # F_na ( S*u ) = S * F_{S^-1*na} ( u )
    for a in range(len(ws)):    # atom index in super cell
        # compute equivalent atom location from inverse rotation in WS cell 
        # temp = S^-1 * a

        # inverse rotate, then to fractional (the other way does not work...)
        # invL * rot * ws[a] == (invL * rot^-1 * L) * (invL * ws[a])
        temp = numpy.dot(invL.T, numpy.dot(numpy.linalg.inv(rot), ws[a]))
        
        # find nb so that b = S^-1 * a: only on the equivalent atoms 'rot'
        for b in range(len(ws)):
            found1 = False
            tmp1  = numpy.dot(invL.T, ws[b]) # to fractional
            # check that the fractional part is the same
            delta = temp - tmp1
            if numpy.linalg.norm(delta - numpy.rint(delta)) < symprec:
                #          ws[b] = rot^-1 * ws[a] : 'b' is equivalent to 'a' by 'rot'
                #    rot * ws[b] =          ws[a]
                #           F[b] = rot^-1 * F[a] 
                #    rot *  F[b] = F[a] 
                found1 = True
                # the rotation matrix rot is applied as is to
                # the force in PHON/set_forces.    
                nforce[a] = numpy.dot(rot, force[b])  # apply rotation
                break # for b
    
        if not found1:
            nforce = None
            break # for a
        # end for a
    
    return nforce

# ------------------------------------------------------------------------------  
def _phonons_run_force2(phonon, dxlist, force1, a, sign, symprec=1e-6):  
    """From a given force set, we derive the forces in Cartesian coordinates
       
       code outrageously adapted from PHON/src/set_forces by D. Alfe
       
       input:
          phonon: ASE Phonon object with Atoms and Calculator
          dxlist: displacement list, in Cartesian coordinates
          force1: the forces list (3 moves) determined for dxlist
          a:      the index of the atom being moved
          sign:   the [+-] directions to move ([0,1] -> [+,-])
          symprec:the precision for comparing positions
          
       output:
          The forces for Cartesian displacements along [xyz]
    
    """

    found = 0
    # print 'dxlist:', dxlist
    # print 'force1:', force1
    # get the 3 displacement vectors (rows) and normalise
    dxnorm   = numpy.asarray([dx/numpy.linalg.norm(dx) for dx in dxlist])
    invdx    = numpy.linalg.inv(dxnorm)
    identity = numpy.diag([1,1,1])
        
    if numpy.linalg.norm(invdx + identity) < symprec: # inv(dx) = -I -> I
        invdx=identity
    # project the 3 forces into that Cartesian basis
    force2 = []
    for i in range(len(dxlist)):
        force2.append(invdx[i,0] * force1[0] \
                    + invdx[i,1] * force1[1] \
                    + invdx[i,2] * force1[2])
    # write the pickle files assuming we have now [xyz]
    for i in range(3):      # xyz
        # skip if the pickle already exists
        filename = '%s.%d%s%s.pckl' % \
               (phonon.name, a, 'xyz'[i], ' +-'[sign])
        if os.path.isfile(filename):
            # Skip if already done. Also the case for initial 'disp/dx'
            continue
        # write the pickle for the current Cartesian axis
        fd = opencew(filename)
        if rank == 0:
            pickle.dump(force2[i], fd, protocol=2)
            sys.stdout.write('Writing %s\n' % filename)
            fd.close()
            sys.stdout.flush()
        found += 1
    
    if found > 0 and False:
        print '[ASE] Displacements, and inverse:'
        print dxlist
        print invdx
        print 'force2 [xyz]:', force2
    return force2
    
# ------------------------------------------------------------------------------
def phonons_read(phonon, method='Frederiksen', symmetrize=3, acoustic=True,
         cutoff=None, born=False, **kwargs):
    """Read forces from pickle files and calculate force constants.

    Extra keyword arguments will be passed to ``read_born_charges``.
    
    This implementation is similar to the ASE one, but can make use of different
    gradient estimates, depending on what is available on disk (pickles).
    Can use:
      displacement .[xyz]+
      displacement .[xyz]-
      equilibrium  .qe
      
    
    Parameters
    ----------
    method: str
        Specify method for evaluating the atomic forces.
    symmetrize: int
        Symmetrize force constants (see doc string at top) when
        ``symmetrize != 0`` (default: 3). Since restoring the acoustic sum
        rule breaks the symmetry, the symmetrization must be repeated a few
        times until the changes a insignificant. The integer gives the
        number of iterations that will be carried out.
    acoustic: bool
        Restore the acoustic sum rule on the force constants.
    cutoff: None or float
        Zero elements in the dynamical matrix between atoms with an
        interatomic distance larger than the cutoff.
    born: bool
        Read in Born effective charge tensor and high-frequency static
        dielelctric tensor from file.
        
    """

    # proceed with pure ASE 'Phonon' object.
    method = method.lower()
    assert method in ['standard', 'frederiksen']
    if cutoff is not None:
        cutoff = float(cutoff)
        
    # Read Born effective charges and optical dielectric tensor
    if born:
        phonon.read_born_charges(**kwargs)
    
    # Number of atoms
    natoms = len(phonon.indices)
    # Number of unit cells
    N = numpy.prod(phonon.N_c)
    # Matrix of force constants as a function of unit cell index in units
    # of eV / Ang**2
    C_xNav = numpy.empty((natoms * 3, N, natoms, 3), dtype=float)
    
    # get equilibrium forces (if any)
    filename = phonon.name + '.eq.pckl'
    feq = 0
    if os.path.isfile(filename):
        feq = pickle.load(open(filename))
        if method == 'frederiksen':
            for i, a in enumerate(phonon.indices):
                feq[a] -= feq.sum(0)

    # Loop over all atomic displacements and calculate force constants
    for i, a in enumerate(phonon.indices):
        for j, v in enumerate('xyz'):
            # Atomic forces for a displacement of atom a in direction v
            basename = '%s.%d%s' % (phonon.name, a, v)
            
            if os.path.isfile(basename + '-.pckl'):
                fminus_av = pickle.load(open(basename + '-.pckl'))
            else:
                fminus_av = None
            if os.path.isfile(basename + '+.pckl'):
                fplus_av = pickle.load(open(basename + '+.pckl'))
            else:
                fplus_av = None
            
            if method == 'frederiksen': # translational invariance
                if fminus_av is not None:
                    fminus_av[a] -= fminus_av.sum(0)
                if fplus_av is not None:
                    fplus_av[a]  -= fplus_av.sum(0)
            
            if fminus_av is not None and fplus_av is not None:
                # Finite central difference derivative
                C_av = (fminus_av - fplus_av)/2
            elif fminus_av is not None:
                # only the - side is available: forward difference
                C_av =  fminus_av - feq
            elif fplus_av is not None:
                # only the + side is available: backward difference
                C_av = -(fplus_av - feq)

            C_av /= phonon.delta  # gradient

            # Slice out included atoms
            C_Nav = C_av.reshape((N, len(phonon.atoms), 3))[:, phonon.indices]
            index = 3*i + j
            C_xNav[index] = C_Nav

    # Make unitcell index the first and reshape
    C_N = C_xNav.swapaxes(0 ,1).reshape((N,) + (3 * natoms, 3 * natoms))

    # Cut off before symmetry and acoustic sum rule are imposed
    if cutoff is not None:
        phonon.apply_cutoff(C_N, cutoff)
        
    # Symmetrize force constants
    if symmetrize:
        for i in range(symmetrize):
            # Symmetrize
            C_N = phonon.symmetrize(C_N)
            # Restore acoustic sum-rule
            if acoustic:
                phonon.acoustic(C_N)
            else:
                break
         
    # Store force constants and dynamical matrix
    phonon.C_N = C_N
    phonon.D_N = C_N.copy()
    
    # Add mass prefactor
    m_a = phonon.atoms.get_masses()
    phonon.m_inv_x = numpy.repeat(m_a[phonon.indices]**-0.5, 3)
    M_inv = numpy.outer(phonon.m_inv_x, phonon.m_inv_x)
    for D in phonon.D_N:
        D *= M_inv

# ------------------------------------------------------------------------------
# compatibility with PhonoPy
# ------------------------------------------------------------------------------
def find_primitive(cell, symprec=1e-5):
    """
    A primitive cell is searched in the input cell. When a primitive
    cell is found, an object of Atoms class of the primitive cell is
    returned. When not, None is returned.
    
    From phonopy/structure/symmetry.py
    """
    
    # return as is when spglib is not installed
    if not has_spglib:
        return cell
        
    lattice, positions, numbers = spglib.find_primitive(cell, symprec)
    if lattice is None:
        return cell
    else:
        return Atoms(numbers=numbers,
                     scaled_positions=positions,
                     cell=lattice,
                     pbc=True)


def phonopy_run(phonon, single=True, filename='FORCE_SETS'):
    """Run the phonon calculations, using PhonoPy.
    
    input:
    
      phonon: ASE Phonon object with Atoms and Calculator
      single: when True, the forces are computed only for a single step, and then
              exit. This allows to split the loop in independent iterations. When
              calling again the 'run' method, already computed steps are ignored,
              missing steps are completed, until no more are needed. When set to
              False, all steps are done in a row.
    
    output:
    
      True when a calculation step was performed, False otherwise or no more is needed.
      the force constants are then stored in the Phonon ASE object.

    """
    
    from phonopy import Phonopy
    from phonopy.structure.atoms import Atoms as PAtoms
    from phonopy.structure.atoms import PhonopyAtoms
    import phonopy.file_IO as file_IO
    
    # we first look if an existing phonon pickle exists. This is the case if we
    # are running with iterative calls while return value is True. The first call
    # will then create the objects, which are subsequently updated until False.
    
    # Set calculator if provided
    assert phonon.calc is not None, "Provide calculator in Phonon __init__ method"
    
    # Atoms in the supercell -- repeated in the lattice vector directions
    # beginning with the last
    supercell = phonon.atoms * phonon.N_c
    
    # create a PhonopyAtoms object
    cell=PhonopyAtoms(phonon.atoms.get_chemical_symbols(), 
        positions=phonon.atoms.get_positions(), 
        cell=phonon.atoms.get_cell(), magmoms=None)

    # is there an existing PhonoPy calculation ?
    # using factor=6.46541380e-2=VaspToeV
    if os.path.exists('FORCE_SETS'):
        phonpy = Phonopy(cell, numpy.diag(phonon.N_c), 
            is_auto_displacements=False,
            primitive_matrix= None,
            dynamical_matrix_decimals= None,
            force_constants_decimals= None,
            symprec= 1e-05,
            is_symmetry= True,
            use_lapack_solver= False,
            log_level= 1)
        force_sets = file_IO.parse_FORCE_SETS(filename='FORCE_SETS')
        phonpy.set_displacement_dataset(force_sets)
        # inactivate magmoms in supercell as some calculators do not provide that
        phonpy._supercell.magmoms=None
        phonpy.produce_force_constants(calculate_full_force_constants=False)
    else:
        # create a PhonoPy Phonon object.
        phonpy = Phonopy(cell, numpy.diag(phonon.N_c))
        # generate displacements (minimal set)
        phonpy.generate_displacements(distance=0.01)
        # iterative call for all displacements
        set_of_forces, flag = phonopy_run_calculate(phonon, phonpy, supercell, single)
        
        if flag is True:
            return flag # some more work is probably required
            
        sys.stdout.write('[ASE/Phonopy] Computing force constants\n')
        # use symmetry to derive forces in equivalent displacements
        phonpy.produce_force_constants(forces=set_of_forces)
        
        # generate disp.yaml and FORCE_SETS (for later use)
        displacements = phonpy.get_displacements()
        directions    = phonpy.get_displacement_directions()
        file_IO.write_disp_yaml(displacements,
                                phonpy.get_supercell(),
                                directions=directions)
        file_IO.write_FORCE_SETS(phonpy.get_displacement_dataset())
    
    # store as additional data in atoms 'info'
    phonon.atoms.info["phonopy"] = phonpy
    
    # save the PhonoPy object
    fid = opencew("phonopy.pkl")
    if fid is not None and rank == 0:
        print "[ASE/Phonopy] Writing %s" % "phonopy.pkl"
        pickle.dump(phonpy, fid, protocol=2)
        fid.close()
    

    # transfer results to the ASE phonon object
    # Number of atoms (primitive cell)
    natoms = len(phonon.indices)
    # Number of unit cells (supercell)
    N = numpy.prod(phonon.N_c)
    
    # Phonopy: force_constants size is [N*natoms,N*natoms,3,3]
    # Phi[i,j,a,b] with [i,j = atom in supercell] and [a,b=xyz]
    force_constants = phonpy.get_force_constants()
    # the atoms [i] which are moved are in the first cell of the supercell, i.e.Ni=0
    # the forces are then stored for all atoms [Nj,j] as [3,3] matrices
    
    # we compute the sum on all supercells, which all contain n atoms.
    C_N = numpy.zeros((N, 3*natoms, 3*natoms), dtype=complex)
    Ni=0
    for Nj in range(N):
        for ni in range(natoms):
            Nni = ni
            for nj in range(natoms):
                # compute Nn indices
                Nnj = Nj*natoms + nj
                # get fc 3x3 matrix
                C_N[Nj,(3*ni):(3*ni+3),(3*nj):(3*nj+3)] += force_constants[Nni][Nnj]
            
    # convert to ASE storage
    # ASE: phonon.C_N size is be [N, 3*natoms, 3*natoms]
    # Phi[i,j] = Phi[j,i]
    phonon.C_N = C_N
    
    # fill dynamical matrix (mass prefactor)
    phonon.D_N = phonon.C_N.copy()
        
    # Add mass prefactor
    m_a = phonon.atoms.get_masses()
    phonon.m_inv_x = numpy.repeat(m_a[phonon.indices]**-0.5, 3)
    M_inv = numpy.outer(phonon.m_inv_x, phonon.m_inv_x)
    for D in phonon.D_N:
        D *= M_inv
    
    return False  # nothing left to do
    
# ------------------------------------------------------------------------------
def phonopy_run_calculate(phonon, phonpy, supercell, single):
    # calculate forces using PhonoPy: phonopy/example/ase/8Si-phonon.py
    
    # get the displacements
    disps = phonpy.get_displacements()
    
    # get displaced supercells: api_phonopy._build_supercells_with_displacements()
    supercells = phonpy.get_supercells_with_displacements()
        
    # Do calculation on equilibrium structure (used to improve gradient accuracy)
    supercell.set_calculator(phonon.calc)
    if rank == 0:
        print "[ASE/Phonopy] Computing equilibrium"
    feq = phonons_run_eq(phonon, supercell)
    for i, a in enumerate(phonon.indices):
        try:
            feq[a] -= feq.sum(0)  # translational invariance
        except AttributeError:
            pass
    
    # compute the forces
    set_of_forces   = []
    for d in range(len(supercells)):
        # skip if this step has been done already.
        filename = '%s.%d.pkl' % (phonon.name, d)
        if os.path.isfile(filename):
            f = open(filename, 'rb')
            forces = pickle.load(f)
            f.close()
        else:
            # proceed with force calculation for current displaced supercell
            disp = disps[d]
            scell = supercells[d]
            if rank == 0:
                print "[ASE/Phonopy] Computing step %i/%i" % (d+1, len(supercells))
                print "Moving  atom #%-3i %-3s    to " % \
                    (disp[0], scell.get_chemical_symbols()[disp[0]]), disp[1:]
            
            cell = Atoms(symbols=scell.get_chemical_symbols(),
                         scaled_positions=scell.get_scaled_positions(),
                         cell=scell.get_cell(),
                         pbc=True)
            cell.set_calculator(phonon.calc)
            forces = cell.get_forces()
            
            drift_force = forces.sum(axis=0)
            # print "[Phonopy] Drift force:", "%11.5f"*3 % tuple(drift_force)
            # Simple translational invariance
            for force in forces:
                force -= drift_force / forces.shape[0]
            
            if feq is not None:
                try:
                    forces -= feq # forward difference, but assumes equilibrium is not always 0
                except AttributeError:
                    print "[ASE/PhonoPy] Can not use forward difference (equilibrium forces are mis-formatted)"

            # save the forces in a pickle
            f = opencew(filename)
            if f is not None and rank == 0:
                sys.stdout.write('Writing %s\n' % filename)
                pickle.dump(forces, f, protocol=2)
                f.close()
            
            # in single shot mode, we return when forces were computed
            if single:
                return set_of_forces, True
        
        # store the incremental list of forces
        set_of_forces.append(forces)
        
    return set_of_forces, False
    
# ------------------------------------------------------------------------------
def phonopy_band_structure(phonpy, path_kc, modes=False):
    """
    Calculate phonon dispersion along a path in the Brillouin zone, using PhonoPy

        The dynamical matrix at arbitrary q-vectors is obtained by Fourier
        transforming the real-space force constants. In case of negative
        eigenvalues (squared frequency), the corresponding negative frequency
        is returned.

        Eigenvalues and modes are in units of PhonoPy and Ang/sqrt(amu),
        respectively. The default PhonoPy energy unit is THz.
        
        This implementation provides better eigenvector estimates that the ASE one.

        Parameters
        ----------
        path_kc: ndarray
            List of k-point coordinates (in units of the reciprocal lattice
            vectors) specifying the path in the Brillouin zone for which the
            dynamical matrix will be calculated.
        modes: bool
            Returns both frequencies and modes when True.
    """
         
    D = phonpy._dynamical_matrix
    num_atom     = len(D._p2s_map)
    # pre-allocating arrays brings a speed improvement
    omega_kl     = numpy.zeros((len(path_kc), 3 * num_atom))
    u_kl         = numpy.zeros((len(path_kc), 3 * num_atom, num_atom, 3), dtype=complex)
    # the call to phonopy._set_dynamical_matrix() is expensive. 
    # It create a  DynamicalMatrix() at each call, but this has been done 
    # in phonopy.produce_force_constants(). So we avoid it.
    
    # pre compute masses sqrt(mass * mass)
    m_inv_x = numpy.repeat(numpy.sqrt(D._pcell.get_masses()), 3)
    mass    = numpy.outer(m_inv_x, m_inv_x)
    
    for iqc, q_c in enumerate(path_kc):
    
        D.set_dynamical_matrix(q_c)
        dm = D.get_dynamical_matrix()
        
        if modes:
            omega2_l, u_xl = numpy.linalg.eigh(dm, UPLO='U')
            # Sort eigenmodes according to eigenvalues (see below) and
            # multiply with mass prefactor
            u_lx = (m_inv_x[:, numpy.newaxis] *
                    u_xl[:, omega2_l.argsort()]).T.copy()
            u_kl[iqc] = u_lx.reshape((-1, num_atom, 3))
        else:
            omega2_l = numpy.linalg.eigvalsh(dm, UPLO='U')
                
        # Sort eigenvalues in increasing order
        omega2_l.sort()
        # Use dtype=complex to handle negative eigenvalues
        omega_l = numpy.sqrt(omega2_l.astype(complex))

        # Take care of imaginary frequencies
        if not numpy.all(omega2_l >= 0.):
            indices = numpy.where(omega2_l < 0)[0]
            
            omega_l[indices] = -1 * numpy.sqrt(numpy.abs(omega2_l[indices].real))

        omega_kl[iqc]=omega_l.real
        
    if modes:
        return numpy.asarray(omega_kl) * phonpy._factor, numpy.asarray(u_kl)
    else:
        return numpy.asarray(omega_kl) * phonpy._factor
        
