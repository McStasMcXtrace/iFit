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
    from ase.spacegroup import Spacegroup
except ImportError:
    from ase.lattice.spacegroup import Spacegroup

import sys
import pickle
import numpy

from ase.atoms import Atoms
from ase.utils import opencew
from ase.parallel import rank, barrier

# ------------------------------------------------------------------------------
def get_spacegroup(atoms, symprec=1e-5):
    """Determine the spacegroup to which belongs the Atoms object.
    
    Parameters:
    atoms:    an Atoms object
    symprec:  Symmetry tolerance, i.e. distance tolerance in Cartesian 
              coordinates to find crystal symmetry.
              
    The Spacegroup object is returned, and stored in atoms.info['spacegroup'] 
    when this key does not exist (avoids overwrite). To force overwrite of the 
    spacegroup first use:
        del atoms.info["spacegroup"]
                  
    Examples:
    
    >>> from ase.lattice import bulk
    >>> atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    >>> sg = get_spacegroup.get_spacegroup(atoms)
    """ 
    #
    # we try all available spacegroups from 1 to 230
    # a Space group is the collection of all symmetry operations which let the 
    # unit cell invariant.
    found    = None
    positions  = atoms.get_scaled_positions(wrap=True)
    positions -= positions[0] # make sure we are insensitive to translation
    
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
def phonons_run(phonon):
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

    """

    # Atoms in the supercell -- repeated in the lattice vector directions
    # beginning with the last
    atoms_N = phonon.atoms * phonon.N_c
    
    # Set calculator if provided
    assert phonon.calc is not None, "Provide calculator in __init__ method"
    atoms_N.set_calculator(phonon.calc)
    
    # Do calculation on equilibrium structure
    filename = phonon.name + '.eq.pckl'

    fd = opencew(filename)
    if fd is not None:
        # Call derived class implementation of __call__
        output = phonon.__call__(atoms_N)
        # Write output to file
        if rank == 0:
            pickle.dump(output, fd)
            sys.stdout.write('Writing %s\n' % filename)
            fd.close()
        sys.stdout.flush()
        # test if the equilibrium is good
        is_equilibrated = numpy.max(output) < 1e-10
        if not is_equilibrated:
            warnings.warn('WARNING: The lattice is not equilibrated ! The force estimate may be wrong.')

    # Positions of atoms to be displaced in the reference cell
    natoms = len(phonon.atoms)
    offset = natoms * phonon.offset
    pos = atoms_N.positions[offset: offset + natoms].copy()
    
    # Loop over all displacements
    for a in phonon.indices:
        for i in range(3):
            for sign in [-1, 1]:
                # Filename for atomic displacement
                filename = '%s.%d%s%s.pckl' % \
                           (phonon.name, a, 'xyz'[i], ' +-'[sign])
                # Wait for ranks before checking for file
                # barrier()
                fd = opencew(filename)
                if fd is None:
                    # Skip if already done
                    continue

                # Update atomic positions
                atoms_N.positions[offset + a, i] = \
                    pos[a, i] + sign * phonon.delta
                
                # Call derived class implementation of __call__
                output = phonon.__call__(atoms_N)
                # Write output to file
                if rank == 0:
                    pickle.dump(output, fd)
                    sys.stdout.write('Writing %s\n' % filename)
                    fd.close()
                sys.stdout.flush()
                
                # fill equivalent displacements from spacegroup
                _phonons_run_symforce(phonon, atoms_N, \
                    atoms_N.positions[offset + a], output, pos)
                
                # Return to initial positions
                atoms_N.positions[offset + a, i] = pos[a, i]
                
# ------------------------------------------------------------------------------
def _phonons_run_symforce(phonon, atoms, disp, force, pos):
    """From a given force set, we derive the forces for equivalent displacements
    by applying the corresponding symmetry operators."""
    
    # check if a spacegroup is defined
    if 'spacegroup' not in atoms.info or atoms.info["spacegroup"] is None:
        return
    
    pos0 = pos.copy() # store the current equilibrium positions
    
    # get the equivalent displacements for the current move
    # the 'sites' are wrapped, and normalized.
    disp,kinds,rot,trans = equivalent_sites(atoms.info["spacegroup"], \
                disp, onduplicates='keep')
    
    # Loop over all planned displacements (past and future)
    # and search for such a move that matches one of the equivalent 'disp'
    for a in phonon.indices:
        for i in range(3):
            for sign in [-1, 1]:
                # skip if the pickle already exists
                filename = '%s.%d%s%s.pckl' % \
                       (phonon.name, a, 'xyz'[i], ' +-'[sign])
                fd = opencew(filename)
                if fd is None:
                    # Skip if already done
                    continue
                
                # compute the 'planned' displacements (atom, xyz, sign)
                pos = pos0.copy()
                pos[a, i] = pos[a, i] + sign * phonon.delta
                # wrap this move, so that we can compare
                pos = equivalent_sites(atoms.info["spacegroup"], pos[a], \
                    onduplicates='keep')[0]
                pos = pos[0]  # one single displaced atom
                
                # does this move belongs to the equivalent displacements ?
                if pos in disp:
                    # if found in sites, we apply rot/trans to forces
                    # get the index of the matching site for rot and trans
                    index = numpy.where((disp == pos).all(axis=1))
                    try:
                        index=index[0][0]
                    except IndexError:
                        raise IndexError('Failed to apply symmetry operators. File %s is left empty.' % filename)
                    # apply rotation/translation on forces
                    nforce = force.copy()
                    for j in range(len(nforce)):
                        nforce[j] = numpy.dot(rot[index], nforce[j]) + trans[index]
                    # write the pickle
                    if rank == 0:
                        pickle.dump(nforce, fd)
                        sys.stdout.write('Writing %s (from spacegroup "%s" symmetry)\n' % (filename, atoms.info["spacegroup"].symbol))
                        fd.close()
                
