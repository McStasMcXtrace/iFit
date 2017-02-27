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
    positions -= positions[0] 
    
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

                # Update atomic positions. Shift is applied to Cartesian coordinates.
                disp    = [0,0,0]
                disp[i] = sign * phonon.delta # in Angs, cartesian
                atoms_N.positions[offset + a, i] = \
                    pos[a, i] + disp[i]
                
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
                    disp, output, a)
                
                # Return to initial positions
                atoms_N.positions[offset + a, i] = pos[a, i]
                
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

