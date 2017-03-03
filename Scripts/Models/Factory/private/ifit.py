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

# from ase.lattice import bulk

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

from ase.atoms import Atoms
from ase.utils import opencew
from ase.parallel import rank, barrier

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
def get_spacegroup(atoms, symprec=1e-5, centre=None):
    """Determine the spacegroup to which belongs the Atoms object.
    
    Parameters:
    atoms:    an Atoms object
    symprec:  Symmetry tolerance, i.e. distance tolerance in Cartesian 
              coordinates to find crystal symmetry.
    centre:   None, or the index of the atoms to use as centre
              
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
    
    # we try all available spacegroups from 1 to 230, backwards
    # a Space group is the collection of all symmetry operations which lets the 
    # unit cell invariant.
    found    = None
    positions  = atoms.get_scaled_positions(wrap=True)  # in the lattice frame
    
    # make sure we are insensitive to translation. this choice is arbitrary and 
    # could lead to a 'slightly' wrong guess for the Space group, e.g. do not  
    # guess centro-symmetry.
    if centre:
        try:
            positions -= positions[centre]
        except IndexError:
            pass
    
    # search space groups from the highest symmetry to the lowest
    # retain the first match
    for nb in range(230,0,-1):
        sg        = Spacegroup(nb)
        #
        # now we scan all atoms in the cell and look for equivalent sites
        sites,kinds,rot,trans = equivalent_sites(sg, positions, 
                onduplicates='keep', symprec=symprec)
        #    
        # the equivalent sites should match all other atom locations in the cell
        # as the spacegroup transforms the unit cell in itself
        if len(sites) == len(positions):
            # store the space group into the list
            found = sg
            break
    #
    # return None when no space group is found (would be surprising)
    if found is not None and 'spacegroup' not in atoms.info:
          atoms.info["spacegroup"] = found
    #   
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
        rots: list
            A list of rotation operators to apply in order to produce equivalent sites
        transs: list
            A list of translation operators to apply in order to produce equivalent sites

        Example:

        >>> from ase.lattice.spacegroup import Spacegroup
        >>> sg = Spacegroup(225)  # fcc
        >>> sites, kinds,rot,trans = sg.equivalent_sites([[0, 0, 0], [0.5, 0.0, 0.0]])
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
        rots  = []
        transs= []
        
        scaled = numpy.array(scaled_positions, ndmin=2)
        for kind, pos in enumerate(scaled):
            for rot, trans in sg.get_symop():
                # check if R.x + T = x' is in the atoms list (equivalent site)
                # x is expressed in fractional/scaled/direct/lattice coordinates
                site = numpy.mod(numpy.dot(rot, pos) + trans, 1.)
                if not sites:
                    sites.append(site)
                    kinds.append(kind)
                    rots.append(rot)
                    transs.append(trans)
                    continue
                t = site - sites
                mask = numpy.all((abs(t) < symprec) |
                              (abs(abs(t) - 1.0) < symprec), axis=1)
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
                    rots.append(rot)
                    transs.append(trans)
        return numpy.array(sites), kinds, rots, transs
        
# ------------------------------------------------------------------------------
def phonons_run(phonon, single=True, usesymmetry=False):
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
    
      phonon object with Atoms and Calculator
    
    output:
    
      True when a calculation step was performed, false otherwise

    """

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
    if usesymmetry:
        signs = [1]
    else:
        signs = [-1,1]
    for a in phonon.indices:
        for i in range(3):
            for sign in signs:
                # Update atomic positions. Shift is applied to Cartesian coordinates.
                disp    = [0,0,0]
                disp[i] = sign * phonon.delta # in Angs, cartesian
                
                # we determine the point group, i.e. we center the atom 'a'
                atoms0 = phonon.atoms.copy()                   # initial lattice, not displaced
                atoms0.positions      -= atoms0.positions[a] # centre on 'a'
                
                sg = get_spacegroup(atoms0)
                
                # we get the scaled positions for the atoms
                atoms0.positions[a,i] += disp[i]
                scaled = atoms0.get_scaled_positions() # this is used by the Calculator
                
                # skip this move if it is equivalent to an existing one using the
                # space group rotations
                # filename, rotation = _phonon_run_equiv(sg, phonon, atoms0.get_cell(), scaled[a])
                # if filename:
                #    continue
                
                # Filename for atomic displacement
                filename = '%s.%d%s%s.pckl' % \
                           (phonon.name, a, 'xyz'[i], ' +-'[sign])
                # Wait for ranks before checking for file
                # barrier()
                fd = opencew(filename)
                if fd is None:
                    # Skip if already done
                    continue
                
                # positions is in Cartesian coordinates, but this is e.g. 
                # converted to fractional coordinates when calling the calculator
                # in ase.calculator: FileIOCalculator.calculate -> write_input
                atoms_N.positions[offset + a, i] = \
                    pos[a, i] + disp[i]
                    
                scaled = atoms_N.get_scaled_positions()
                print "Moving atom #%i %s by " % (a, atoms_N.get_chemical_symbols()[a]), disp, " (Angs)\n"
                print "Fractional coordinates ", scaled[a],"\n"
                
                # Call derived class implementation of __call__
                output = phonon.__call__(atoms_N)
                # Write output to file
                if rank == 0:
                    pickle.dump(output, fd)
                    sys.stdout.write('Writing %s\n' % filename)
                    fd.close()
                sys.stdout.flush()
                
                # fill equivalent displacements from spacegroup
                if usesymmetry and False:
                    _phonons_run_symforce(phonon, atoms_N, \
                       disp, output, a)
                
                # Return to initial positions
                atoms_N.positions[offset + a, i] = pos[a, i]
                
                if single:
                    return True # and some more iterations may be required

    return False  # nothing left to do

# ------------------------------------------------------------------------------       
def _phonon_run_equiv(sg, phonon, cell, move):
    """Check if the displacement 'disp' matches an existing move when applying
       the space group.
       
       Parameters:
       
       sg:    space group 
       cell:  atoms.get_cell
       disp:  scaled (fractional) displacement in 'atoms'
       
       Returns:
       
       filename: the filename pickle that holds the equivalent move
       rotation: the rotation matrix to convert filename/forces into 'disp'
    """
    
    sites, kinds, rotations, translations = equivalent_sites(sg, move)
    
    L    = cell                 # lattice cell            = at
    invL = numpy.linalg.inv(L)  # no 2*pi multiplier here = inv(at) = bg'
    
    # now scan all planned displacements
    for a in phonon.indices:
        for i in range(3):
            for sign in [-1, 1]:
            
                # the 'planned' displacement
                disp    = [0,0,0]
                disp[i] = sign * phonon.delta # in Angs, cartesian
                
                # convert from Cartesian to scaled positions
                scaled = numpy.mod(numpy.dot(invL, disp), 1.)
                
                # check if the corresponding file exists
                # Filename for atomic displacement
                filename = '%s.%d%s%s.pckl' % \
                           (phonon.name, a, 'xyz'[i], ' +-'[sign])
                
                # is it in the equivalent 'sites' and filename exists ?
                if scaled in sites and os.path.isfile(filename):
                    # get the index of the matching equivalent move
                    index = numpy.where(numpy.isclose(sites,scaled).all(axis=1))
                    if len(index):
                        index = index[0]
                    if len(index) > 1:
                        print "WARNING: ifit._phonon_run_equiv: multiple equivalent sites ! Using first out of ", index
                    if len(index):                    
                        index=index[0]
                    else:
                        continue  # no symmetry operator found
    
                    return filename, rotations[index]
    return None, None
    
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
                # Finite difference derivative
                C_av = (fminus_av - fplus_av)/2
            elif fminus_av is not None:
                # only the - side is available
                C_av =  fminus_av
            elif fplus_av is not None:
                # only the + side is available
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
                
# ------------------------------------------------------------------------------
def _phonons_run_symforce(phonon, atoms, disp, force, a):
    """From a given force set, we derive the forces for equivalent displacements
    by applying the corresponding symmetry operators."""
    
    # check if a spacegroup is defined
    if 'spacegroup' not in atoms.info or atoms.info["spacegroup"] is None:
        return
        
    L    = atoms.get_cell()     # lattice cell            = at
    invL = numpy.linalg.inv(L)  # no 2*pi multiplier here = inv(at) = bg'

    # get the equivalent displacements for the current move
    # the first row is: rot=I and trans=0
    
    # equivalent_sites uses the scaled positions, so we must convert the 'disp'
    # from Cartesian to lattice frame
    disp_lat = numpy.dot(invL, disp)
    disp_lat = numpy.mod(disp_lat, 1.)
    
    # search for all [disp_lat*rot+trans] equivalents
    disps_lat,kinds,rot,trans = equivalent_sites(atoms.info["spacegroup"], \
                disp_lat, onduplicates='keep')
    
    # **** Coordinate frames ***************************************************
    # fractional/lattice coordinate use the [a,b,c] frame, and coordinates are 
    #   fractions in the cell. The cell frame may not be ortho-normal. This frame
    #   is also labeled as 'Direct'.
    #   atoms.scaled_positions is using the Direct/lattice frame
    #       fractional = np.linalg.solve(self.cell.T, self.positions.T).T
    #                  = inv(self.cell) * positions
    # Cartesian coordinates use an ortho-normal frame.
    #   atoms.positions is using the cartesian frame
    #       positions = np.dot(scaled_positions, self._cell)
    #   displacements in ASE are in the cartesian frame
    
    # **** Transformations *****************************************************
    # rot[index] applied on disp gives new_disp equivalent site, when given
    # in lattice frame => rot[index] is then the 'lattice' rotation matrix.
    # Then we conclude on rotation operators:
    #   lattice:   rot[index]         e.g. a simple matrix
    #   cartesian: R= numpy.dot(L, numpy.dot(rot[index], invL))
    # displacements:
    #   cartesian: cart   = L*direct, e.g. [-.01 0 0] a 'nice' vector xyz in ASE
    #   direct:    direct = invL*cart e.g. a mixed vector

    # Loop over all planned displacements (past and future)
    # and search for such a move that matches one of the equivalent 'disp'
    for i in range(3):
        for sign in [-1, 1]:
            
            # compute the 'planned' displacements (atom, xyz, sign)
            new_disp    = [0,0,0]
            new_disp[i] = sign * phonon.delta # in Angs, cartesian
            # convert this displacement from cartesian to lattice frame
            new_disp_lat= numpy.dot(invL, new_disp)
            # modulo in the cell, Direct/scaled/lattice/fractional
            new_disp_lat= numpy.mod(new_disp_lat, 1.)
            
            # does this move belongs to the equivalent displacements ?
            if new_disp_lat in disps_lat:
                # if found in sites, we apply rotation to forces
                # get the index of the matching site for rot and trans
                index = numpy.where(numpy.isclose(disps_lat,new_disp_lat).all(axis=1))
                if len(index):
                    index = index[0]
                if len(index) > 1:
                    print "WARNING: ifit:phonons:run_symforce: multiple equivalent sites ! Using first out of ", index
                if len(index):                    
                    index=index[0]
                else:
                    continue  # no symmetry operator found
                    
                # skip if the pickle already exists
                filename = '%s.%d%s%s.pckl' % \
                       (phonon.name, a, 'xyz'[i], ' +-'[sign])
                fd = opencew(filename)
                if fd is None:
                    # Skip if already done. Also the case for initial 'disp'
                    continue

                nforce = force.copy()

                # convert rotation operator from lattice to cartesian frame
                # L x R x L^-1 as in phononpy.similarity_transformation
                R      = numpy.dot(L, numpy.dot(rot[index], invL))
                
                # apply rotation operator in cartesian frame on forces
                nforce = numpy.dot(R, nforce.T).T

                # write the pickle
                if rank == 0:
                    pickle.dump(nforce, fd)
                    sys.stdout.write('Writing %s (from spacegroup "%s" symmetry)\n' % (filename, atoms.info["spacegroup"].symbol))
                    fd.close()

