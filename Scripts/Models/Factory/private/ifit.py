""" This module contains functions to be use for Phonon calculations.
  
  get_spacegroup:
    define an equivalent to spglib.get_spacegroup
    Try out with:
        from ase.lattice import bulk
        import get_spacegroup
        atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
        sg = get_spacegroup.get_spacegroup(atoms)
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

# ------------------------------------------------------------------------------
def get_spacegroup(atoms, symprec=1e-5):
    """Determine the spacegroup to which belongs the Atoms object.
    
    Parameters:
    atoms:    an Atoms object
    symprec:  Symmetry tolerance, i.e. distance tolerance in Cartesian 
              coordinates to find crystal symmetry.
    center:   None, or the index of the atoms to use as center
              
    The Spacegroup object is returned, and stored in atoms.info['spacegroup'] 
    when this key does not exist (avoids overwrite). To force overwrite of the 
    spacegroup first use:
        del atoms.info["spacegroup"]
                  
    Examples:
    
    >>> from ase.lattice import bulk
    >>> atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    >>> sg = get_spacegroup.get_spacegroup(atoms)
    """

    # use spglib when it is available (and return)
    if has_spglib:
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
        sites,kinds = equivalent_sites(sg, positions, 
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
def equivalent_sites(sg, scaled_positions, onduplicates='keep',
                         symprec=1e-3):
        """Returns the scaled positions and all their equivalent sites.

        Parameters:

        scaled_positions: list | array
            List of non-equivalent sites given in unit cell coordinates.
        onduplicates : 'keep' | 'replace' | 'warn' | 'error'
            Action if `scaled_positions` contain symmetry-equivalent
            positions:
            
            'keep'
               ignore additional symmetry-equivalent positions
            'replace'
                replace
            'warn'
                like 'keep', but issue an UserWarning
            'error'
                raises a SpacegroupValueError
                    
        symprec: float
            Minimum "distance" betweed two sites in scaled coordinates
            before they are counted as the same site.

        Returns:

        sites: array
            A NumPy array of equivalent sites.
        kinds: list
            A list of integer indices specifying which input site is
            equivalent to the corresponding returned site.

        Example:

        >>> from ase.lattice.spacegroup import Spacegroup
        >>> sg = Spacegroup(225)  # fcc
        >>> sites, kinds = sg.equivalent_sites([[0, 0, 0], [0.5, 0.0, 0.0]])
        >>> sites
        array([[ 0. ,  0. ,  0. ],
               [ 0. ,  0.5,  0.5],
               [ 0.5,  0. ,  0.5],
               [ 0.5,  0.5,  0. ],
               [ 0.5,  0. ,  0. ],
               [ 0. ,  0.5,  0. ],
               [ 0. ,  0. ,  0.5],
               [ 0.5,  0.5,  0.5]])
        >>> kinds
        [0, 0, 0, 0, 1, 1, 1, 1]
        """
        kinds = []
        sites = []
        
        scaled = numpy.array(scaled_positions, ndmin=2)
        for kind, pos in enumerate(scaled):
            for rot, trans in sg.get_symop():
                # check if R.x + T = x' is in the atoms list (equivalent site)
                # x is expressed in fractional/scaled/direct/lattice coordinates
                
                site = numpy.mod(numpy.dot(rot, pos) + trans, 1.)
                if not sites:
                    sites.append(site)
                    kinds.append(kind)
                    continue
                t = site - sites
                # test if this is a new site location, and append it
                mask = numpy.linalg.norm(t - numpy.rint(t) < symprec)
                if numpy.any(mask):
                    ind = numpy.argwhere(mask)[0][0]
                    if kinds[ind] == kind:
                        pass
                    elif onduplicates == 'keep':
                        pass
                    elif onduplicates == 'replace':
                        kinds[ind] = kind
                    elif onduplicates == 'warn':
                        warnings.warn('scaled_positions %d and %d '
                                      'are equivalent'%(kinds[ind], kind))
                    elif onduplicates == 'error':
                        raise SpacegroupValueError(
                            'scaled_positions %d and %d are equivalent'%(
                                kinds[ind], kind))
                    else:
                        raise SpacegroupValueError(
                            'Argument "onduplicates" must be one of: '
                            '"keep", "replace", "warn" or "error".')
                else:
                    sites.append(site)
                    kinds.append(kind)
        return numpy.array(sites), kinds
        
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
                
                if force0 is not None and rot is not None:
                    # we will now apply the rotation to the force array
                    output = _phonons_run_force1_rot(force0, ws, invL, rot, symprec=1e-6)
                    if output is None:
                        # failed using symmetry to derive force. Trigger full computation.
                        force0 = None
                    else:
                        print "Imaging atom #%-3i %-3s    to " % \
                            (offset + a, supercell.get_chemical_symbols()[a]), pos[a] + disp, \
                            " (Angs) using rotation:"
                    print rot
                if force0 is None or rot is None: # force0 is None or rot is None
                    # move atom 'a' by 'disp'
                    supercell.positions[offset + a] = pos[a] + disp
                    print "Moving  atom #%-3i %-3s    to " % \
                        (offset + a, supercell.get_chemical_symbols()[a]), pos[a] + disp, " (Angs)"
                        
                    # Call derived class implementation of __call__
                    output = phonon.__call__(supercell)
                    
                    # Return to initial positions
                    supercell.positions[offset + a] = pos[a]
                    
                    # store the forces for subsequent moves using symmetry
                    force0 = output
                
                # append the forces to the force1 list
                if output is None:
                    print "Warning: force1 is None !!"
                    
                print "ASE/Phonon: forces are: ", output.shape

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
       The test is for colinearity, and singularity of the formed basis
       
       input:
           dxlist:  a list of vectors
           dx:      a new vector which is tested
           symprec: precision for comparing vectors
       output:
           True when dx is independent from dxlist, False otherwise.
    """

    # test if the rotated displacement is collinear to
    # a stored one (in list 'dxlist'). test is done on nomalised vectors.
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
            pickle.dump(force2[i], fd)
            sys.stdout.write('Writing %s\n' % filename)
            fd.close()
            sys.stdout.flush()
        found += 1
    
    if found > 0 and False:
        print 'Displacements, and inverse:'
        print dxlist
        print invdx
        print 'force2 [xyz]:', force2
    return force2
    
# ------------------------------------------------------------------------------
def phonon_read(phonon, method='Frederiksen', symmetrize=3, acoustic=True,
         cutoff=None, born=False, **kwargs):
    """Read forces from pickle files and calculate force constants.

    Extra keyword arguments will be passed to ``read_born_charges``.
    
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
            
            if method == 'frederiksen':
                if fminus_av is not None:
                    fminus_av[a] -= fminus_av.sum(0)
                if fplus_av is not None:
                    fplus_av[a]  -= fplus_av.sum(0)
            
            if fminus_av is not None and fplus_av is not None:
                # Finite central difference derivative
                C_av = (fminus_av - fplus_av)/2
            elif fminus_av is not None:
                # only the - side is available: forward difference
                C_av =  fminus_av
            elif fplus_av is not None:
                # only the + side is available: backward difference
                C_av = -fplus_av

            C_av /= phonon.delta  # gradient

            # Slice out included atoms
            C_Nav = C_av.reshape((N, len(phonon.atoms), 3))[:, phonon.indices]
            index = 3*i + j
            C_xNav[index] = C_Nav

    # Make unitcell index the first and reshape
    print "ASE/Phonon: force constants C_xNav are: ", C_xNav.shape
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
    print "ASE/Phonon: force constants are: ", C_N.shape
    
    # Add mass prefactor
    m_a = phonon.atoms.get_masses()
    phonon.m_inv_x = numpy.repeat(m_a[phonon.indices]**-0.5, 3)
    M_inv = numpy.outer(phonon.m_inv_x, phonon.m_inv_x)
    for D in phonon.D_N:
        D *= M_inv
                

# ------------------------------------------------------------------------------
# compatibility with PhonoPy
# ------------------------------------------------------------------------------

def phonon_run_phonopy(phonon, single=True):
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
    from phonopy.structure.atoms import Atoms as PhonopyAtoms
    
    # we first look if an existing phonon pickle exists. This is the case if we
    # are running with iterative calls while return value is True. The first call
    # will then create the objects, which are subsequently updated until False.
    
        
    # Atoms in the supercell -- repeated in the lattice vector directions
    # beginning with the last
    supercell = phonon.atoms * phonon.N_c
    
    # Set calculator if provided
    assert phonon.calc is not None, "Provide calculator in Phonon __init__ method"
    
    # create a PhonopyAtoms object
    bulk=PhonopyAtoms(supercell.get_chemical_symbols(), 
        positions=supercell.get_positions(), 
        cell=supercell.get_cell())

    # create a PhonoPy Phonon object. The supercell is already in the 'bulk'
    phonpy = Phonopy(bulk,[[1,0,0],[0,1,0],[0,0,1]])
    
    # generate displacements (minimal set)
    phonpy.generate_displacements(distance=0.01)
    disps = phonpy.get_displacements()
    for d in disps:
        print "[ASE/Phonopy]", d[0], d[1:]
    
    # get displaced supercells
    supercells = phonpy.get_supercells_with_displacements()
    
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
            scell = supercells[d]
            cell = Atoms(symbols=scell.get_chemical_symbols(),
                         scaled_positions=scell.get_scaled_positions(),
                         cell=scell.get_cell(),
                         pbc=True)
            cell.set_calculator(phonon.calc)
            forces = cell.get_forces()
            print "ASE/PhonoPy: forces are: ", forces.shape
            
            drift_force = forces.sum(axis=0)
            # print "[Phonopy] Drift force:", "%11.5f"*3 % tuple(drift_force)
            # Simple translational invariance
            for force in forces:
                force -= drift_force / forces.shape[0]

            # save the forces in a pickle
            f = open(filename, 'wb')
            pickle.dump(forces, f)
            f.close()
            
            # in single shot mode, we return when forces were computed
            if single:
                return True
        
        # store the incremental list of forces
        set_of_forces.append(forces)
        
    # use symmetry to derive forces in equivalent displacements
    phonpy.produce_force_constants(forces=set_of_forces)
    
    # transfer results to the ASE phonon object
    # Number of atoms (primitive cell)
    natoms = len(phonon.indices)
    # Number of unit cells (supercell)
    N = numpy.prod(phonon.N_c)
    
    # Phonopy: force_constants size is [N*natoms,N*natoms,3,3]
    # which is Phi(i,j,a,b) with [i,j = atom in supercell] and [a,b=xyz]
    force_constants = phonpy.get_force_constants()
    print "ASE/PhonoPy: force constants are: ", force_constants.shape, natoms
    force_constants = numpy.reshape(force_constants, (N,N,3 * natoms, 3 * natoms))
    print "ASE/PhonoPy: force constants are now: ", force_constants.shape, (N,N,3 * natoms, 3 * natoms)
    
    # ASE: phonon.C_N size is be [N, 3*natoms, 3*natoms]
    phonon.C_N = force_constants[0,:,:,:]
    print "ASE/PhonoPy: C_N is: ", phonon.C_N.shape
    
    # fill dynamical matrix (mass prefactor)
    phonon.D_N = phonon.C_N.copy()
        
    # Add mass prefactor
    m_a = phonon.atoms.get_masses()
    phonon.m_inv_x = numpy.repeat(m_a[phonon.indices]**-0.5, 3)
    M_inv = numpy.outer(phonon.m_inv_x, phonon.m_inv_x)
    for D in phonon.D_N:
        D *= M_inv
    
    return False  # nothing left to do
