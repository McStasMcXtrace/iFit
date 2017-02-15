function ret = sqw_phon2ase(target)
  % get a directory containing a PHON/PhonoPy output, and create necessary items so that
  % it can be used by ASE
  
  ret = []; flag_phon = true; flag_phonopy = true;
  % look if this is a PHON output directory: INPHON, POSCAR
  if isempty(dir(fullfile(target,'FORCES'))) || isempty(dir(fullfile(target,'DISP')))
    flag_phon = false;
  end
  if isempty(dir(fullfile(target,'FORCE_SET')))
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
