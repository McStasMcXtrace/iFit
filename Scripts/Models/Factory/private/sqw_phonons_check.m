function [options, result, read] = sqw_phonons_check(configuration, options, status)
% sqw_phonons_check: read the initial material structure
%   requires: nothing except file/directory with material 'POSCAR' or other
%   output:   generates a script 'sqw_phonons_check.py' to generate ASE atoms.pkl

target = options.target;

% determine if the atoms.pkl exists. If so, nothing else to do
if ~isempty(dir(fullfile(target, 'atoms.pkl')))
  disp([ mfilename ': re-using ' fullfile(target, 'atoms.pkl') ]);
  return
end

if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
else           precmd = ''; end

if isdir(configuration)
  % search for 'known' configurations in the directory from PHON/PhonoPy
  files = search_files(target, { 'POSCAR*','SPOSCAR' });
  flag_get_supercell = true;
  supercell = [];
  if ~isempty(files)
    files = files.name;
    if strncmp(files, 'POSCAR',6)
      configuration = fullfile(target, files);
    elseif strncmp(files, 'SPOSCAR',7)
      configuration = fullfile(target, files);
      supercell     = [ 1 1 1 ];
      flag_get_supercell = false;
    end
  end
else flag_get_supercell=false;
end

if flag_get_supercell
  files = search_files(target, ...
    { 'phonon.yaml','INPHON','quasiharmonic_phonon.yaml','band.yaml'});
  if ~isempty(files)
    files = iLoad(fullfile(target, files));
  end
  if isfield(files.Data, 'NDIM')
    supercell = files.Data.NDIM;
  elseif isfield(files.Data, 'supercell_matrix')
    supercell = trace(files.Data.supercell_matrix);
  end
  if ~isempty(supercell), options.supercell = supercell; 
  elseif any(options.supercell == 0)
    disp([ mfilename ': ERROR: unspecified supercell for previous computation.' ]); 
  end
end

% handle input configuration: read
if exist(configuration)
  read = sprintf('import ase.io\nconfiguration = "%s"\natoms = ase.io.read(configuration)\n', ...
    configuration);
elseif ischar(configuration)
  read = configuration;
  % ASE has changed some of the modules hierarchy from 3.9 to 3.10+
  switch strtok(configuration, ' (')
  case 'bulk'
    % ASE 3.9 and 3.10+
    read = sprintf('from ase.lattice import bulk\natoms = %s\n', configuration);
  case 'molecule'
    if status.ase <= 3.9
      % ASE 3.9
      read = sprintf('from ase.structure import molecule\natoms = %s\n', configuration);
    else
      % ASE 3.10+: from ase.build import molecule
      read = sprintf('from ase.build import molecule\natoms = %s\n', configuration);
    end
  case 'nanotube'
    if status.ase <= 3.9
      % ASE 3.9
      read = sprintf('from ase.structure import nanotube\natoms = %s\n', configuration);
    else
      % ASE 3.10+: from ase.build import nanotube
      read = sprintf('from ase.build import nanotube\natoms = %s\n', configuration);
    end
  case 'crystal'
    if status.ase <= 3.9
      % ASE 3.9
      read = sprintf('from ase.lattice.spacegroup import crystal\natoms = %s\n', configuration);
    else
      % ASE 3.10+: from ase.spacegroup import crystal
      read = sprintf('from ase.spacegroup import crystal\natoms = %s\n', configuration);
    end
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
  'try:    sg = ifit.get_spacegroup(atoms)\n' ...
  'except: sg = None\n' ...
  'fid = open("atoms.pkl","wb")\n' ...
  'pickle.dump(atoms, fid)\n' ...
  'fid.close()' ...
  '# export atoms into usual formats\n' ...
  'print "Exporting structure to usual formats...\\n"\n' ...
  'from ase.io import write\n' ...
  'atomsN = atoms*tuple([2,2,2]) # nicer rendering for viewers\n' ...
  'try: write("configuration.png", atomsN)\n' ...
  'except: pass\n', ...
  'write("configuration.eps", atomsN)\n' ...
  'write("configuration.pov", atomsN)\n' ...
  'write("configuration.cif", atoms, "cif")\n' ...
  'write("configuration.x3d", atomsN, "x3d")\n' ...
  'write("configuration.pdb", atoms, "pdb")\n' ...
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
  [st, result] = system([ precmd 'python ' fullfile(target,'sqw_phonons_check.py') ]);
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
