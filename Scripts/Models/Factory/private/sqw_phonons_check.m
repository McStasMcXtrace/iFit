function [options, result, read] = sqw_phonons_check(configuration, options, status)
% sqw_phonons_check: read the initial material structure
%   requires: nothing except file/directory with material 'POSCAR' or other

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
  'fid = open("' fullfile(target, 'atoms.pkl') '","wb")\n' ...
  'pickle.dump(atoms, fid)\n' ...
  'fid.close()' ...
  '# export atoms into usual formats\n' ...
  'print "Exporting structure to usual formats...\\n"\n' ...
  'from ase.io import write\n' ...
  'write("' fullfile(target, 'configuration.png') '", atoms)\n' ...
  'write("' fullfile(target, 'configuration.eps') '", atoms)\n' ...
  'write("' fullfile(target, 'configuration.pov') '", atoms)\n' ...
  'write("' fullfile(target, 'configuration.cif') '", atoms, "cif")\n' ...
  'write("' fullfile(target, 'configuration.x3d') '", atoms, "x3d")\n' ...
  'write("' fullfile(target, 'configuration.pdb') '", atoms, "pdb")\n' ...
  'write("' fullfile(target, 'configuration.html') '", atoms, "html")\n' ...
  'write("' fullfile(target, 'configuration.etsf') '", atoms, "etsf")\n' ...
  'write("' fullfile(target, 'configuration_SHELX.res') '", atoms, "res")\n' ...
  'write("' fullfile(target, 'configuration_POSCAR') '", atoms, "vasp")\n' ... 
  'import scipy.io as sio\n' ...
  'print "Exporting structure properties...\\n"\n' ...
  'properties = {\n' ...
  '  "reciprocal_cell" : atoms.get_reciprocal_cell()*6.283185, \n' ...
  '  "cell"            : atoms.get_cell(), \n' ...
  '  "volume"          : atoms.get_volume(), \n' ...
  '  "chemical_formula": atoms.get_chemical_formula(), \n' ...
  '  "chemical_symbols": atoms.get_chemical_symbols(), \n' ...
  '  "cell_scaled_unit": atoms.get_scaled_positions(), \n' ...
  '  "masses"          : atoms.get_masses(), \n' ...
  '  "positions"       : atoms.get_positions(), \n' ...
  '  "atomic_numbers"  : atoms.get_atomic_numbers() }\n' ...
  'try:\n' ...
  '  import spglib\n' ...
  'except ImportError:\n' ...
  '  try:\n' ...
  '    from pyspglib import spglib\n' ...
  '  except ImportError:\n' ...
  '    pass\n' ...
  'try:\n' ...
  '  properties["spacegroup"] = s    = spglib.get_spacegroup(atoms)\n' ...
  '  properties["spacegroup_number"] = int(s[s.find("(")+1:s.find(")")])\n' ...
  'except:\n' ...
  '  properties["spacegroup_number"] = properties["spacegroup"] = None\n' ...
  '# export properties as pickle\n' ...
  'fid = open("' fullfile(target, 'properties.pkl') '","wb")\n' ...
  'pickle.dump(properties, fid)\n' ...
  'fid.close()\n' ...
  '# export properties as MAT\n' ...
  'sio.savemat("' fullfile(target, 'properties.mat') '", properties)\n' ...
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

% display message at start
disp(' ')
disp([ mfilename ': starting phonons computation [' datestr(now) ']' ])
if ~isdeployed && usejava('jvm') && usejava('desktop')
  disp([ '  directory  = <a href="' target '">' target '</a>' ]);
else
  disp([ '  directory  = ' target ]);
end
disp(sprintf([ '  material   = ' configuration ]));
disp(sprintf([ '  calculator = ' options.calculator ]));
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
