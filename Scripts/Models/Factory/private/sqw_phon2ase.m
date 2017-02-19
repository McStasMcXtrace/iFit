function ret = sqw_phon2ase(target)
  % get a directory containing a PHON/PhonoPy output, and create necessary items so that
  % it can be used by ASE
  
  % we have to reconstruct the Phonon.get_force_constant=Phonon.C_N set in Phonon.read
  
  ret = []; flag_phon = true; flag_phonopy = true;
  % look if this is a PHON output directory: INPHON, POSCAR
  if isempty(dir(fullfile(target,'FORCES'))) || isempty(dir(fullfile(target,'DISP')))
    flag_phon = false;
  elseif isempty(dir(fullfile(target,'FORCE_SET')))
    flag_phonopy = false;
  end
  if ~flag_phon && ~flag_phonopy
    return
  end
  % look for the initial material POSCAR
  flag_get_supercell = false;
  if ~isempty(dir(fullfile(target,'phonon.pkl')))
    % model is fully prepared by ASE. nothing to do.
    ret = 'phonon.pkl';
    return
  elseif ~isempty(dir(fullfile(target,'atoms.mat')))
    % model is prepared by ASE. Use that.
    ret = 'atoms.mat';
    flag_get_supercell = true;
  elseif ~isempty(dir(fullfile(target,'SPOSCAR')))
    % model not available. Reconstruct it.
    % if SPOSCAR and use supercell=1
    [options, result, read] = sqw_phonons_check(fullfile(target,'SPOSCAR'), options, status);
    supercell               = [ 1 1 1 ];
  elseif ~isempty(dir(fullfile(target,'POSCAR')))
    % if POSCAR and look in INPHON for NDIM=supercell
    [options, result, read] = sqw_phonons_check(fullfile(target,'POSCAR'), options, status);
    flag_get_supercell = true;
  elseif ~isempty(dir(fullfile(target,'POSCAR-unitcell')))
    % if POSCAR-unitcell
    [options, result, read] = sqw_phonons_check(fullfile(target,'POSCAR-unitcell'), options, status);
    flag_get_supercell = true;
  end
  % look for the supercell
  if flag_get_supercell
    % phonon.yaml from PhonoPy: supercell_matrix
    % INPHON      from PHON:    NDIM
  end
  
  % look for the displacements
  % disp.yaml from PhonoPy
  % DISP      from PHON
  
% check 
  
% * create a Phonons object, with qe_ase calculator when available, or None
if status.quantumespresso_ase
  [decl, calc, signal] = sqw_phonons_calc(options, status, 'quantumespresso_ase', read);
else
  calc = 'None';  % used in ph.run(), so not needed
end

% * ph=Phonons(atoms, supercell=(n,n,n)); 
% sprintf('ph = Phonons(atoms, calc, supercell=(%i, %i, %i), delta=0.05)
% * compute what is missing from the read FORCES
%    'C_N = FORCES; '
%    '# Store force constants and dynamical matrix'
%    'ph.C_N = C_N'
%    'ph.D_N = C_N.copy()'
%    '# Add mass prefactor'
%    'm_a = ph.atoms.get_masses()'
%    'ph.m_inv_x = np.repeat(m_a[ph.indices]**-0.5, 3)'
%    'M_inv = np.outer(ph.m_inv_x, ph.m_inv_x)'
%    'for D in ph.D_N:'
%    '  D *= M_inv'
% * save it into phonon.pkl'
%    'fid = open("' fullfile(target, 'phonon.pkl') '","wb")\n' , ...
%    'calc = ph.calc\n', ...
%    sav, ...
%    'pickle.dump(ph, fid)\n', ...
%    'fid.close()\n', ...
% * the python eval script calls, Phonons.band_structure and Phonons.dos
% * requies ph.D_N, ph.N_c = supercell, 

% ------------------------------------------------------------------------------

function [FORCES, natoms] = sqw_phon_forceset(lines)
% convert a FORCE_SET/PhonoPy file into a FORCES/PHON file

% remove all empty lines
index = cellfun('isempty', lines);
lines(index) = [];

% read the number of atoms in the supercell
natoms = str2double(lines{1});
if ~isscalar(natoms), forces = []; return; end

% read the number of displacements
ndisp  = str2num(lines{2});
if ~isscalar(ndisp), forces = []; return; end

displacements= zeros(ndisp, 4);
forces       = cell(ndisp,1);

% skip empty lines and read
next_line_disp = 0; disp_index=0; next_line_force=0; force_index=0;
for index=3:numel(lines)
  nb = str2num(lines{index}); % current line
  
  if next_line_force  % when reading the forces after a displacement vector
    force_cat(next_line_force, :) = nb;
    next_line_force = next_line_force +1;
    if next_line_force==natoms,
      next_line_force=0; % stop reading this displacement forces
      forces{force_index} = force_cat;
    end
  end
  
  if next_line_disp % when reading a displacement vector, then activate forces read
    disp_index     = disp_index + 1;
    next_line_disp = 0;
    displacements(disp_index, :) = [ atom_moved nb ];
    next_line_force= 1;
  end
  
  if isscalar(nb) % single value: index of atom which is moved
    atom_moved = nb; % the atom id which is moved
    % next line should contain the displacement vector
    next_line_disp = 1;
    force_index = force_index+1;
    force_cat = zeros(natoms, 3);
  end
  
end

% now write the FORCES as a string
  FORCES = sprintf('%i\n', numel(forces));  % nb of displacements
  for index=1:numel(forces)
    FORCES = [ FORCES sprintf('%s\n', num2str(displacements(index,:))) ];
    this = cellstr(num2str(forces{index}));
    FORCES = [ FORCES sprintf('    %s\n', this{:}) ];
  end

% ------------------------------------------------------------------------------

% PHONOPY generates FORCE_SET, but only for reduced displacements
% N -> nb of atoms in supercell
% D -> nb of displacements = nb of blocks below
%   ID of atom being moved
%   DISP -> displacement vector xyz
%   FORCES -> (at) rows of (xyz) force components in supercell, e.g. 
%   for Al and 333 supercell (27 atoms), this is a 27x3 matrix.
% PHONOPY generates POSCAR-unitcell and phonon.yaml
% YAML.read('phonon.yaml') -> supercell_patrix
% ase.read(POSCAR) -> atoms


function read_phon_forces(target)
% PHON generates DISP as a reduced set of vectors (taking into account symmetries)
% FORCES also contains the displacement vectors
% N -> nb of displacements/block (taking into account symmetries)
% (at) (displacement vector xyz)
% (at) rows of (xyz) force components in supercell, e.g. 
% for Al and 333 supercell (27 atoms), this is a 27x3 matrix.
% ------------------------------------------------------------------------------

function read_ase_forces(target)
% ASE generates pickle files phonon.(at)(xyz)(+-).pckl
% each pickle file has (at) rows of (xyz) force components, e.g. 
% for Al and 333 supercell (27 atoms), this is a 27x3 matrix.

