function [options, result, read] = sqw_phonons_check(configuration, options, status)
% sqw_phonons_check: read the initial material structure
%   requires: nothing except file/directory with material 'POSCAR' or other
%   output:   generates a script 'sqw_phonons_check.py' to generate ASE atoms.pkl

target = options.target;
result = [];
read = [];

% determine if the atoms.pkl exists. If so, nothing else to do
if ~isempty(dir(fullfile(target, 'atoms.pkl')))
  disp([ mfilename ': re-using ' fullfile(target, 'atoms.pkl') ]);
end

if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
else           precmd = ''; end

config_dir         = [];
% handle import of full directory: search POSCAR, forces, past pickles, ...
if isdir(configuration)
  config_dir = configuration;
  % search for 'known' configurations in the given directory from PHON/PhonoPy
  file = search_files(configuration, { 'atoms.pkl', 'POSCAR*','SPOSCAR*','*_POSCAR'});
  
  if ~isempty(file)
    copyfile(fullfile(config_dir, file.name), target);
    disp([ mfilename ': Re-using lattice cell ' file.name ' from ' config_dir ]);
    configuration = fullfile(target, file.name);
    if strncmp(file, 'SPOSCAR',7)
      options.supercell     = [ 1 1 1 ];
    end
  end
end

if ~isempty(config_dir)
  % get any previous supercell definition
  file = search_files(config_dir, ...
    { 'phonopy.yaml','phonon.yaml','INPHON','quasiharmonic_phonon.yaml','band.yaml'});

  if ~isempty(file) && any(options.supercell == 0)
    file = iLoad(fullfile(config_dir, file.name));
    
    if ~isempty(file)
        if isfield(file.Data, 'NDIM')
          supercell = file.Data.NDIM;
        elseif isfield(file.Data, 'supercell_matrix')
          supercell = diag(file.Data.supercell_matrix);
        else supercell = [];
        end
        if ~isempty(supercell)
          options.supercell = supercell;
          disp([ mfilename ': Re-using supercell ' mat2str(options.supercell) ]); 
        elseif any(options.supercell == 0)
          disp([ mfilename ': WARNING: unspecified supercell from previous computation.' ]); 
        end
    end
  end

  % get FORCES from previous PhonoPy calculation
  for f={'FORCE_SETS','disp.yaml'}
    % look if there is are previous PhonoPy files, and copy them to the target
    if exist(fullfile(config_dir, f{1}))
      try
        copyfile(fullfile(config_dir, f{1}), target);
        disp([ mfilename ': Re-using ' f{1} ' from ' config_dir ]); 
      end
    end
  end
  
  % get previous ASE phonon pickle files
  files = dir(fullfile(config_dir, 'phonon.*.*'));
  for f=files(:)'
    copyfile(fullfile(config_dir, f.name), target);
    disp([ mfilename ': Re-using ' f.name ' from ' config_dir ]); 
  end
  % get previous ASE vibration pickle files
  files = dir(fullfile(config_dir, 'vib.*.*'));
  for f=files(:)'
    copyfile(fullfile(config_dir, f.name), target);
    disp([ mfilename ': Re-using ' f.name ' from ' config_dir ]); 
  end
end



% handle input configuration: read
if exist(configuration)
  [~,~,e] = fileparts(configuration);
  if strcmp(e, '.pkl')
    read = sprintf('import pickle\nconfiguration = "%s"\natoms = pickle.load(open(configuration,"rb"))\n', configuration);
  else
    read = sprintf('import ase.io\nconfiguration = "%s"\natoms = ase.io.read(configuration)\n', ...
      configuration);
  end
elseif ischar(configuration)
  read = configuration;
  % ASE has changed some of the modules hierarchy from 3.9 to 3.10+
  switch strtok(configuration, ' (')
  case 'bulk'
    read = sprintf([ 'try:\n' ...
     '    from ase.build import bulk\n' ...
     'except ImportError:\n' ...
     '    from ase.lattice import bulk\n' ...
     'atoms = %s\n' ], configuration);
  case 'molecule'
    read = sprintf([ 'try:\n' ...
      '    from ase.build import molecule\n' ...
      'except ImportError:\n' ...
      '    from ase.structure import molecule\n' ...
      'atoms = %s\n' ], configuration);
  case 'nanotube'
    read = sprintf([ 'try:\n' ...
      '    from ase.build import nanotube\n' ...
      'except ImportError:\n' ...
      '    from ase.structure import nanotube\n' ...
      'atoms = %s\n' ], configuration);
  case 'crystal'
    read = sprintf([ 'try:\n' ...
      '    from ase.spacegroup import crystal\n' ...
      'except ImportError:\n' ...
      '    from ase.lattice.spacegroup import crystal\n' ...
      'atoms = %s\n' ], configuration);
  end
end

options.script_create_atoms = sprintf([ ...
  '# python script built by ifit.mccode.org/Models.html sqw_phonons\n' ...
  '# on ' datestr(now) '\n' ...
  '# E. Farhi, Y. Debab and P. Willendrup, J. Neut. Res., 17 (2013) 5\n' ...
  '# S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002.\n' ...
  '#\n' ...
  '# read initial material structure and save atoms as a pickle\n\n' ...
  read ...
  'import pickle\n' ...
  'from os import chdir\n' ...
  'import ifit\n' ...
  'chdir("' target '")\n' ...
  'try:\n' ...
  '    sg = ifit.get_spacegroup(atoms)\n' ...
  '    atoms1 = ifit.find_primitive(atoms)\n' ...
  '    if atoms1:\n' ...
  '        atoms = atoms1\n' ...
  'except: sg = None\n' ...
  'fid = open("atoms.pkl","wb")\n' ...
  'pickle.dump(atoms, fid)\n' ...
  'fid.close()' ...
  '# export atoms into usual formats\n' ...
  'print "Exporting structure to usual formats...\\n"\n' ...
  'from ase.io import write\n' ...
  'atomsN = atoms*tuple([2,2,2]) # nicer rendering for viewers\n' ...
  'try: write("configuration.png", atomsN)\n' ...
  'except (ValueError,TypeError): pass\n', ...
  'write("configuration.eps", atomsN)\n' ...
  'write("configuration.pov", atomsN)\n' ...
  'write("configuration.cif", atoms, "cif")\n' ...
  'write("configuration.x3d", atomsN, "x3d")\n' ...
  'try: write("configuration.pdb", atoms, "pdb")\n' ...
  'except ValueError: pass\n', ...
  'write("configuration.html", atomsN, "html")\n' ...
  'try: write("configuration.etsf", atoms, "etsf")\n' ...
  'except: pass\n' ...
  'write("configuration_SHELX.res", atoms, "res")\n' ...
  'write("configuration_POSCAR", atoms, "vasp")\n' ... 
  'import scipy.io as sio\n' ...
  'print "Exporting structure properties...\\n"\n' ...
  'properties = {\n' ...
  '  "reciprocal_cell" : atoms.get_reciprocal_cell()*6.283185307, \n' ...
  '  "cell"            : atoms.get_cell(), \n' ...
  '  "volume"          : atoms.get_volume(), \n' ...
  '  "chemical_formula": atoms.get_chemical_formula(), \n' ...
  '  "chemical_symbols": atoms.get_chemical_symbols(), \n' ...
  '  "cell_scaled_unit": atoms.get_scaled_positions(), \n' ...
  '  "masses"          : atoms.get_masses(), \n' ...
  '  "positions"       : atoms.get_positions(), \n' ...
  '  "atomic_numbers"  : atoms.get_atomic_numbers() }\n' ...
  'if sg is not None:\n' ...
  '  properties["spacegroup"]        = sg.symbol\n' ...
  '  properties["spacegroup_number"] = sg.no\n' ...
  '  properties["spacegroup_rotations"] = sg.get_rotations()\n' ...
  'else:\n' ...
  '  properties["spacegroup_number"] = properties["spacegroup"] = None\n' ...
  '# remove None values in properties\n' ...
  'properties = {k: v for k, v in properties.items() if v is not None}\n' ...
  '# export properties as pickle\n' ...
  'fid = open("properties.pkl","wb")\n' ...
  'pickle.dump(properties, fid)\n' ...
  'fid.close()\n' ...
  '# export properties as MAT\n' ...
  'sio.savemat("properties.mat", properties)\n' ...
  ]);

% we create a python script to read/check/export the initial structure
fid = fopen(fullfile(target,'sqw_phonons_check.py'),'w+');
if fid == -1
  result = 'ERROR fopen';
  sqw_phonons_error([ mfilename ': failed create material structure reader for ' ...
    configuration ' (sqw_phonons_check.py)' ], options);
  return
end
fprintf(fid, '%s\n', options.script_create_atoms);
fclose(fid);

% this script uses our ASE 'extensions' for spacegroup and phonons.run
% we copy this ifit.py as well
if isempty(dir(fullfile(target, 'ifit.py')))
  copyfile(fullfile(fileparts(which(mfilename)),'ifit.py'), target);
end

% default to available core (read from Java)
if ~isfield(options,'mpi') && usejava('jvm')
  r=java.lang.Runtime.getRuntime;
  % mem_avail   = r.freeMemory;
  options.mpi = r.availableProcessors;
end

if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
  options.mpirun = [ status.mpirun ' -n ' num2str(options.mpi) ];
  if isfield(options, 'machinefile')
    options.mpirun = [ options.mpirun ' -machinefile ' options.machinefile ];
  end
end

options.configuration = configuration;

if options.use_phonopy && ...
    (strcmpi(options.calculator,'quantumespresso') || strcmpi(options.calculator,'quantumespresso_phon'))
  options.calculator = 'quantumespresso_ase'; % PhonoPy can NOT be used with PHON
elseif isempty(status.phon) && ...
  (strcmpi(options.calculator,'quantumespresso') || strcmpi(options.calculator,'quantumespresso_phon')) ...
  && ~isempty(status.quantumespresso_ase)
  options.calculator = 'quantumespresso_ase'; % PHON not available -> use QEutil
elseif options.use_phonopy && strcmpi(options.calculator,'gpaw')
  % GPAW can not be used with PhonoPy (issue with write/read PhonoPy pickle)
  options.use_phonopy = 0;
end


% display message at start
disp(' ')
disp([ mfilename ': starting phonons computation [' datestr(now) ']' ])
if ~isdeployed && usejava('jvm') && usejava('desktop')
  disp([ '  directory  = <a href="' target '">' target '</a>' ]);
else
  disp([ '  directory  = ' target ]);
end
disp(sprintf([ '  material   = ' configuration ]));
if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
  disp([ '  calculator = ' options.calculator ' using ' num2str(options.mpi) ' cores' ]);
else
  disp([ '  calculator = ' options.calculator ]);
end
% copy the configuration into the target
if ~isempty(dir(configuration))
  try
  copyfile(configuration, target);
  end
end

% ------------------------------------------------------------------------------
% call python: read/export initial configuration
result = '';
try
  [st, result] = system([ precmd status.python ' ' fullfile(target,'sqw_phonons_check.py') ]);
end
% we test if the pickle file could be written. This way even if the export/save 
% properties fail, we can proceed.
disp(result)
if isempty(dir(fullfile(target, 'atoms.pkl')))  % FATAL
  result = 'ERROR python';
  return
end

% determine optimal kpoints and supercell if left in auto mode
if ~isempty(fullfile(target, 'properties.mat')) && ...
    (any(options.supercell <= 0) || any(options.kpoints <= 0))
  properties = load(fullfile(target, 'properties.mat'));
  if isfield(properties, 'atomic_numbers')
    % get the nb of atoms in the model
    nb_at = numel(properties.atomic_numbers);
    % auto mesh should be 1000/nb_at = k^3*supercell^3 (6750/nb_at for high accuracy)
    if all(options.kpoints > 0)
      supercell = floor((1000/nb_at./prod(options.kpoints))^(1/3));
    else
      supercell = floor((1000/nb_at)^(1/6));
    end
    if any(options.supercell <= 0)
      options.supercell=[ supercell supercell supercell ];
      disp([ '  auto: supercell=' num2str(supercell) ])
    end
    kpoints   = ceil((1000/nb_at./prod(options.supercell))^(1/3));
    if any(options.kpoints <= 0)
      options.kpoints = [ kpoints kpoints kpoints ];
      disp([ '  auto: kpoints=  ' num2str(kpoints) ])
    end
  end
end

% compute the K-points density
if ~isempty(fullfile(target, 'properties.mat'))
  properties = load(fullfile(target, 'properties.mat'));
  if isfield(properties, 'atomic_numbers')
    nb_at = numel(properties.atomic_numbers);
    kpts_density = nb_at.*prod(options.supercell).*prod(options.kpoints);
    disp([ '  kpoints_density=  ' num2str(kpts_density) ])
    if kpts_density < 500
      disp('WARNING: The Monkhorst-Pack grid k-points density is small. Expect low computational accuracy.');
      disp('         Be cautious with results (may show unstable modes/negative/imaginary frequencies.');
    end
  end
end

% ==============================================================================
function files = search_files(target, list)
  % search for a file in the list
  files = '';
  for tosearch = list
    files = dir(fullfile(target, tosearch{1}));
    if isempty(files), continue; end
    files = files(~[files.isdir ]); % skip sub-directories
    if isempty(files), continue; end
    files = files(1); break
  end
