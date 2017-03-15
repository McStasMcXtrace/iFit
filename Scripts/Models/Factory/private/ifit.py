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
    """Our own 'private' implementation of get_spacegroup. This one is less secure
       that the spglib one.
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
def phonons_run(phonon, single=True, usesymmetry=False, difference='central'):
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
      usesymmetry: when True, the symmetry operators are used to minimize the
              amount of force calculations needed.
      difference: the method to use for the force difference (gradient) computation.
              can be 'central','forward','backward'. The central difference is 
              more accurate but requires twice as many force calculations.
    
    output:
    
      True when a calculation step was performed, False otherwise or no more is needed.

    """
    
    # prepare to use symmetry operations
    if usesymmetry:
        # check if a spacegroup is defined, else find it
        if 'spacegroup' not in phonon.atoms.info or phonon.atoms.info["spacegroup"] is None:
            sg    = get_spacegroup(phonon.atoms)

    if difference == 'backward':
        signs = [ 1]     # only compute one side, backward difference
    elif difference == 'forward':
        signs = [-1]     # only compute one side, forward difference
    else:
        signs = [-1,1]   # use central difference

    # Atoms in the supercell -- repeated in the lattice vector directions
    # beginning with the last
    atoms_N = phonon.atoms * phonon.N_c
    
    # Set calculator if provided
    assert phonon.calc is not None, "Provide calculator in __init__ method"
    atoms_N.set_calculator(phonon.calc)

    # Positions of atoms to be displaced in the reference cell
    natoms = len(phonon.atoms)
    offset = natoms * phonon.offset
    pos    = atoms_N.positions[offset: offset + natoms].copy()
    
    # Loop over all displacements  
    for a in phonon.indices:    # atom index to move in cell
    
        if usesymmetry:
            # we determine the Wigner-Seitz cell around atom 'a' (Cartesian)
            ws = _get_wigner_seitz(atoms_N, a, max(phonon.N_c)) # max(phonon.N_c)
        else:
            ws = None
            
        # loop to move atom 'a' in xyz direction and +- side
        for i in range(3):      # xyz
            for sign in signs:  # +-
                # Update atomic positions. Shift is applied to Cartesian coordinates.
                disp     = pos[a].copy()
                disp[i] += sign * phonon.delta # in Angs, Cartesian
                
                # Filename for atomic displacement
                filename = '%s.%d%s%s.pckl' % \
                           (phonon.name, a, 'xyz'[i], ' +-'[sign])
                # Wait for ranks before checking for file
                # barrier()
                fd = opencew(filename)
                if fd is None:
                    # Skip if already done
                    continue
                
                # 'positions' are in Cartesian coordinates, but they are e.g. 
                # converted to fractional coordinates when calling the calculator
                # in ase.calculator: FileIOCalculator.calculate -> write_input
                atoms_N.positions[offset + a, i] = disp[i]
                    
                print "Moving atom #%-3i %-3s    at " % (offset + a, atoms_N.get_chemical_symbols()[a]), disp, " (Angs)"
                
                # scaled = atoms_N.get_scaled_positions() 
                #        = numpy.mod(numpy.dot(invL.T, disp),1.) with L=atoms_N.get_cell()
                # cart   = numpy.dot(atoms_N.get_cell().T, scaled)
                
                # Call derived class implementation of __call__
                output = phonon.__call__(atoms_N)
                # Write output to file
                if rank == 0:
                    pickle.dump(output, fd)
                    sys.stdout.write('Writing %s\n' % filename)
                    fd.close()
                sys.stdout.flush()
                
                # fill equivalent displacements from spacegroup
                if usesymmetry and ws is not None:
                    _phonons_run_symforce(phonon, atoms_N, disp, output, pos, ws, signs)
                
                # Return to initial positions
                atoms_N.positions[offset + a, i] = pos[a, i]
                
                if single:
                    return True # and some more iterations may be required

    return False  # nothing left to do

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
    print "Computing Wigner-Seitz cell..."
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
def _phonons_run_symforce(phonon, atoms, disp, force, pos, ws, signs, symprec=1e-6):
    """From a given force set, we derive the forces for equivalent displacements
       by applying the corresponding symmetry operators.
       
       code outrageously adapted from PHON/src/set_forces by D. Alfe
       
       input:
          phonon: ASE Phonon object with Atoms and Calculator
          atoms:  ASE Atoms object for the supercell
          disp:   displacement vector, in Cartesian coordinates
          scaled: displacement vector, in fractional coordinates
          force:  the forces determined for 'disp'
          pos:    the equilibrium  positions, supercell (Cartesian)
          ws:     the Wigner-Seitz positions, supercell (Cartesian)
          
       output:
          True when a new force was generated from symmetry, False otherwise
    
    """
    
    # check if a spacegroup is defined
    if 'spacegroup' not in atoms.info or atoms.info["spacegroup"] is None:
        return
    
    sg   = atoms.info["spacegroup"] # spacegroup for e.g. spglib
    L    = atoms.get_cell()         # lattice cell            = at
    invL = numpy.linalg.inv(L)      # no 2*pi multiplier here = inv(at) = bg.T
    
    scaled = numpy.mod(numpy.dot(invL.T, disp),1.)      # in fractional coordinates
    
    # scaled = atoms_N.get_scaled_positions() 
    #        = numpy.mod(numpy.dot(invL.T, disp),1.)
    # cart   = numpy.dot(L.T, scaled)
    
    # now we search for a new displacement, and will apply rotation to force
    
    # we try all symmetry operations in the spacegroup
    for rot, trans in sg.get_symop():
    
        # find the equivalent displacement from initial displacement (fractional) and rotation
        dx = numpy.mod(numpy.dot(rot, scaled),1.)     # in fractional: dx = rot * disp

        # the new displacement 'dx' must be one of the planned ones, else continue
        found    = False
        filename = None
        for a in phonon.indices:    # atom index to move in cell
            for i in range(3):      # xyz
                for sign in signs: # +-
                    # check if we already found a rotated displacement
                    if found:
                        break
                        
                    # skip if the pickle already exists
                    filename = '%s.%d%s%s.pckl' % \
                           (phonon.name, a, 'xyz'[i], ' +-'[sign])
                    if os.path.isfile(filename):
                        # Skip if already done. Also the case for initial 'disp/dx'
                        continue
                    
                    # compute the 'planned' displacement (atom, xyz, sign)
                    new_disp     = pos[a].copy()
                    new_disp[i] += sign * phonon.delta # in Angs, Cartesian
                    # convert this displacement from Cartesian to lattice frame
                    new_disp     = numpy.mod(numpy.dot(invL.T, new_disp),1.)
                    delta        = new_disp - dx
                    
                    # compare fractional coordinates: we need the rotated disp ?
                    if numpy.linalg.norm(delta - numpy.rint(delta)) < symprec: 
                        found = True  # and filename holds the missing step
        
        # try an other symmetry operation when this one did not work out
        if not found:
            continue   
    
        # print out new displacement
        print 'Equivalent displacement dx=', numpy.dot(L.T, dx - numpy.rint(dx)), ' (Angs)'
        
        # we will now apply the rotation to the force array
        nforce = force.copy()
    
        # scan all atoms to compute the rotated force matrix
        # F_na ( S*u ) = S * F_{S^-1*na} ( u )
        for a in range(len(ws)):    # atom index to move in cell
            # compute equivalent atom location from inverse rotation in WS cell 
            # temp = S^-1 * a
            
            # temp=numpy.dot(bg.T, numpy.dot(inv(rot),ws[a]))
            # temp = ws[a]                                  # WS coordinate (Cartesian)    
            # temp = numpy.mod(numpy.dot(invL.T, temp),1.)  # to fractional      
            # temp = numpy.dot(rot.T, temp)                 # apply inverse rotation
            temp = numpy.dot(invL.T, numpy.dot(numpy.linalg.inv(rot),ws[a]))
            
            # find nb so that b = S^-1 * a: only on the equivalent atoms 'rot'
            for b in range(len(ws)):
                found1 = False
                # tmp1  = numpy.mod(numpy.dot(invL.T, ws[b]), 1.) # to fractional
                tmp1  = numpy.dot(invL.T, ws[b])
                # check that the fractional part is the same
                delta = temp - tmp1
                if numpy.linalg.norm(delta - numpy.rint(delta)) < symprec:
                    #          ws[b] = rot^-1 * ws[a] : 'b' is equivalent to 'a' by 'rot'
                    #    rot * ws[b] =          ws[a]
                    #           F[b] = rot^-1 * F[a] 
                    #    rot *  F[b] = F[a] 
                    found1 = True
                    # surprisingly, the rotation matrix rot is applied as is to
                    # the force in PHON/set_forces.
                    tmp2   = force[b]
                    #tmp2   = numpy.dot(invL.T, force[b])  # to fractional  
                    tmp2   = numpy.dot(rot, tmp2)         # apply rotation
                    #tmp2   = numpy.dot(L.T, tmp2)         # back to Cartesian
                    nforce[a] = tmp2
                    break # for b
        
            if not found1:
                print "Warning: could not apply symmetry to force [_phonons_run_symforce] for ", filename
        # end for a
                
        # write the pickle for the current 'rot'
        fd = opencew(filename)
        if rank == 0:
            pickle.dump(nforce, fd)
            sys.stdout.write('Writing %s (from spacegroup "%s" symmetry)\n' % (filename, sg.symbol))
            fd.close()
        # and loop to the next symmetry operator for an other displacement

    # end for symop
    return found

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
                
