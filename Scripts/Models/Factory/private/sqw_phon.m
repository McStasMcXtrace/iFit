function signal=sqw_phon(poscar, options)
% model=sqw_phon(poscar, options)  
%
%   iFunc/sqw_phon: computes phonon dispersions from Ab-Initio with QuantumEspresso.
%   A model which computes phonon dispersions from the forces acting between
%     atoms. The input arguments is a VASP-type POSCAR file providing the
%     geometry of the cell (lattice and atom positions). 
%     It is recommended to specify if the system is a metal or an insulator.
%     Additional arguments can be given to control the procedure used to
%     compute a supercell and the forces using QuantumEspresso (ab-initio).
%
%   The simplest usage is to call: sqw_phon('gui') which displays an entry dialog box
%   and proceeds with a fully automatic computation, and plots final results.
%
%   When performing a model evaluation, the DOS is also computed on the same HKL
%     grid and stored in model.UserData.DOS as an iData object. For the evaluation
%     to be correct, the HKL grid should thus be broad, e.g. from 0 to 0.5 
%     on QH,QK and QL, with a rather fine gridding (see example below).
%
% The arguments can be any set of:
% POSCAR: file name to an existing POSCAR
%   a VASP compliant file providing initial system structure, with lattice
%     and position of atoms in a stable, equilibrated configuration.
%   The first line MUST start with the chemical formula of the system, as a 
%     single word, which can be followed by any other comment.
%   The specified file name can be different than POSCAR.
%   If you wish to specify the initial system structure with an other file format
%     and have installed ASE, you should use sqw_phonons(file, ...)
%
% options: a structure with optional settings:
%   options.target =path                   where to store all files and FORCES
%     a temporary directory is created when not specified.
%   options.supercell=scalar or [nx ny nz] supercell size
%   options.disp   =[dx dy dz]             initial displacement direction
%     the initial displacement can e.g. be [ 1 -1 1 ]
%   options.potentials={ cell of strings } list of UPF potentials to use
%     if not given a PBE-PAW potential will be used (when available)
%     these can also be directories, which content is then added to the list.
%   options.kpoints=scalar or [nx ny nz]   Monkhorst-Pack grid
%   options.mpi    =scalar                 number of CPUs to use for PWSCF
%     this option requires MPI to be installed (e.g. openmpi).
%   options.potential_auto=0 or 1          automatic setting of the best PP
%   options.phon   =string                 path to PHON executable, e.g. 'phon'
%     if not given, 'phon' is assumed to be available from PATH.
%   options.command=string                 path to QW/PWSCF executable, e.g. 'pw.x' or 'pw.exe'
%     if not given, 'pw' is assumed to be available from PATH.
%
% options affecting memory usage:
%   options.diagonalization='david' or 'cg'
%     The Davidson method is faster than the conjugate-gradient but uses more
%     memory.
%   options.mixing_ndim=scalar             number of iterations used in mixing
%     default=8. If you are tight with memory, you may reduce it to 4. A larger
%     value will improve the SCF convergence, and use more memory.
%
% options affecting SCF convergence:
%   options.mixing_beta=scalar             mixing factor for self-consistency
%     default=0.7. use 0.3 to improve convergence
%   options.ecutrho=scalar                 kinetic energy cutoff (Ry) for charge
%     density and potential. Default=4*ecutwfc, suitable for PAW. Larger value
%     improves convergence, especially for ultra-soft PP (use 8-12*ecutwfc).
%   options.occupations='metal'            for metals ('smearing') help converge
%                       'insulator'        for insulators
%   options.ecutwfc=scalar                 kinetic energy cutoff (Ry) for
%     wavefunctions. default is 15*natoms in [Ry]. Larger value improves convergence.
%   options.electron_maxstep=scalar        max number of iterations for SCF.
%     default=100. Larger value improves convergence.
%   options.toldfe=scalar                Convergence threshold for 
%     selfconsistency. default=1e-6.
%
% The options can also be entered as a single string with 'field=value; ...'.
%
% Once the model has been created, its use requires that axes are given on
% regular qx,qy,qz grids.
%     
% Example:
%   s=sqw_phon([ ifitpath 'Data/POSCAR_Al'],'metal','mpi=4');
%   qh=linspace(0.01,.5,50);qk=qh; ql=qh; w=linspace(0.01,50,51);
%   f=iData(s,[],qh,qk,ql,w); scatter3(log(f(1,:, :,:)),'filled');
%   figure; plot(s.UserData.DOS); % plot the DOS
%
% WARNING: Single intensity and line width parameters are used here.
%   This model is only suitable to compute phonon dispersions for e.g solid-
%   state materials.
%   PHON(executable 'phon') and QuantumEspresso(executable 'pw') must be installed.
%   The temporary directories (UserData.dir) are not removed.
%
% References: https://en.wikipedia.org/wiki/Phonon
%   PHON: D. Alfe, Computer Physics Communications 180,2622-2633 (2009) 
%     <http://chianti.geol.ucl.ac.uk/~dario>. BSD license.
%   Quantum Espresso: P. Giannozzi, et al J.Phys.:Condens.Matter, 21, 395502 (2009)
%     needed, as the PWSCF (pw.x)  program is used to estimate the FORCES.
%     <http://www.quantum-espresso.org/>. GPL2 license.
%
% input:  p: sqw_phon model parameters (double)
%             p(1)=Amplitude
%             p(2)=Gamma   dispersion DHO half-width in energy [meV]
%             p(3)=Background (constant)
%             p(4)=Temperature of the material [K]
%             p(5)=Energy_scaling. All frequencies are multiplied by this value.
%          or p='guess'
%         qh: axis along QH in rlu (row,double)
%         qk: axis along QK in rlu (column,double)
%         ql: axis along QL in rlu (page,double)
%         w:  axis along energy in meV (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_cubic_monoatomic, sqw_phonons, sqw_sine3d, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>

% Private functions (inline)
%       [poscar, options]=sqw_phon_argin(varargin)
%       [forces, geom] = sqw_phon_forces(geom, options)
%       force = sqw_phon_forces_pwscf(displaced, options)
%       [potentials, potentials_full] = sqw_phon_forces_pwscf_potentials(displaced, options)
%       options = sqw_phon_potentials(geom, options)
%       [phon, pwscf]  = sqw_phon_requirements
%       geom1 = sqw_phon_supercell(poscar, options)
%       forces = sqw_phon_forceset(forces)
%
% Private function (in private)
%       import_poscar
%       export_poscar
%       parse_formula
%       molweight

signal = [];

% ********* handle input arguments *********

[poscar, options] = sqw_phon_argin(poscar, options);
p = fileparts(poscar);
if isempty(p), p=pwd; end

if ~isfield(options,'phon') || isempty(options.phon)
  options.phon = 'phon';
  if ispc, options.phon=[ options.phon '.exe' ]; end
end
if ~isfield(options,'command') || isempty(options.command)
  options.command = 'pw.x';
end

options.configuration = poscar;

% check if we use a pre-computed POSCAR (supercell)/FORCES
if isfield(options,'force_file')
  % FORCES format:
  % 1st line: scalar (=<nb_displacements>)
  % 2nd line: vector4
  % FORCE_SET format:
  % 1st line: scalar (=nb atoms in super cell)
  % 2nd line: scalar (=<nb_displacements>)
  forces = fileread(options.force_file);
  lines  = textscan(forces,'%s','delimiter','\n','whitespace',''); % read all lines
  lines=lines{1}; l1 = str2num(lines{1}); l2 = str2num(lines{2}); 
  % identify if this is a FORCE_SET (PhonoPy) and convert it to FORCES for PHON
  % the number of atoms(1st line) in FORCE_SET must match size(geom.coords,2)
  if isscalar(l1) && isscalar(l2)
    % supposingly FORCE_SET: we convert it to FORCES/PHON
    forces = sqw_phon_forceset(lines);
  end
  % supposingly FORCES, we just store it as is
  signal.UserData.FORCES = forces;
end

% ********* compute the forces *********

% check/create supercell
if ishandle(options.gui), waitbar(0.05, options.gui, [ mfilename ': computing supercell' ]); end
geom = sqw_phon_supercell(poscar, options); %  supercell
if isempty(geom) return; end

% ********* look for potentials *********
if ~isfield(options, 'force_file')
  if ishandle(options.gui), waitbar(0.10, options.gui, [ mfilename ': fetching potentials' ]); end
  
  options = sqw_phon_potentials(geom, options);

  % compute and write FORCES
  if ishandle(options.gui), waitbar(0.15, options.gui, [ mfilename ': compute DynMatrix PHON/QuantumEspresso (be patient)' ]); end
  [forces, geom] = sqw_phon_forces(geom, options);
  if isempty(forces) return; end
  
  signal.UserData.FORCES    = fileread(fullfile(p,'FORCES'));
  signal.UserData.SUPERCELL = fileread(fullfile(p,'POSCAR'));  % supercell POSCAR
else
  geom.potentials = cell(numel(geom.symbols),1);
end

mass = [];
compound = ''; potentials = '';
for index=1:numel(geom.symbols)
  mass(index) = molweight(geom.symbols{index});
  compound    = [ compound geom.symbols{index} num2str(geom.atomcount(index)) ];
  potentials  = [ potentials geom.potentials{index} ' ' ];
end
mass = num2str(mass);

if ishandle(options.gui), waitbar(0.7, options.gui, [ mfilename ': building model' ]); end


% ************** BUILD THE MODEL **************
signal.Name           = [ 'Sqw_phon_' compound ' PHON/QuantumEspresso DHO [' mfilename ']' ];

signal.Description    = [ compound 'S(q,w) 3D dispersion PHON/QuantumEspresso with DHO line shape. Pseudo-Potentials: ' potentials '. Configuration: ' poscar ];

signal.Parameters     = {  ...
  'Amplitude' ...
  'Gamma Damped Harmonic Oscillator width in energy [meV]' ...
  'Background' ...
  'Temperature [K]' ...
  'Energy_scaling' ...
   };
  
signal.Dimension      = 4;         % dimensionality of input space (axes) and result

signal.Guess = [ 1 .1 0 10 1 ];

signal.UserData.POSCAR    = fileread(poscar); % the initial POSCAR file (not supercell)
signal.UserData.dir       = p;
signal.UserData.geometry  = geom;
signal.UserData.options   = options;
signal.UserData.potentials= potentials;
signal.UserData.DOS       = [];

% required to avoid Matlab to use its own libraries
if ismac,      precmd = 'DYLD_LIBRARY_PATH= ; DISPLAY= ; ';
elseif isunix, precmd = 'LD_LIBRARY_PATH= ; DISPLAY= ; '; 
else           precmd=''; end

% for the evauation, in order to use ASE, we need to:
% * create a Phonons object, with qe_ase calculator when available, or None
% * ph=Phonons(atoms, supercell=(n,n,n)); ph.C_N = FORCES; ph.D_N = C_N / M_inv
% * save it into phonon.pkl'
% * the python eval script calls, Phonons.band_structure and Phonons.dos
% * requies ph.D_N, ph.N_c = supercell, 

% get code to read xyzt and build HKL list and convolve DHO line shapes
[script_hkl, script_dho] = sqw_phonons_templates;

if strcmp(options.accuracy,'fast')
  lcentral = 'fprintf(fid, ''  LCENTRAL = .F.\n'');';
else
  lcentral = '';
end

signal.Expression     = { ...
  '% check if FORCES and supercell POSCAR are here', ...
[ 'pw = pwd; target = this.UserData.dir;' ], ...
  'if ~isdir(target), target = tempname; mkdir(target); this.UserData.dir=target; end', ...
[ 'if isempty(dir(fullfile(target, ''SUPERCELL'')))' ], ...
  '  fid=fopen(fullfile(target, ''POSCAR''), ''w'');', ...
  '  if fid==-1, error([ ''model '' this.Name '' '' this.Tag '' could not write POSCAR into '' target ]); end', ...
  '  fprintf(fid, ''%s\n'', this.UserData.SUPERCELL);', ...
  '  fclose(fid);', ...
  'end', ...
[ 'if isempty(dir(fullfile(target, ''FORCES'')))' ], ...
  '  fid=fopen(fullfile(target, ''FORCES''), ''w'');', ...
  '  if fid==-1, error([ ''model '' this.Name '' '' this.Tag '' could not write FORCES into '' target]); end', ...
  '  fprintf(fid, ''%s\n'', this.UserData.FORCES);', ...
  '  fclose(fid);', ...
  'end', ...
  script_hkl{:}, ...
  '% write INPHON for FREQ generation ------------------------------------------', ...
[ 'fid=fopen(fullfile(target,''INPHON''),''w'');' ], ...
[ 'fprintf(fid, ''# Control file for PHON to compute the modes from POSCAR and FORCES in %s\n'', target);' ], ...
  'fprintf(fid, ''#   PHON: <http://chianti.geol.ucl.ac.uk/~dario>\n'');', ...
  'fprintf(fid, ''#   D. Alfe, Computer Physics Communications 180,2622-2633 (2009)\n'');', ...
  'fprintf(fid, ''  LSYMM=.TRUE.\n'');', ...
  'fprintf(fid, ''  LSUPER=.F.\n'');', ...
  lcentral, ...
  'fprintf(fid, ''  LEIGEN=.T.\n'');', ...
[ 'fprintf(fid, ''# NDIM  = ' num2str(options.supercell) '\n'');' ], ...
  'fprintf(fid, ''# LFREE=.T.\n'');', ...
  'fprintf(fid, ''  LFORCEOUT=.TRUE.\n'');', ...
  'fprintf(fid, ''  IPRINT=2\n'');', ...
  'fprintf(fid, ''# TEMPERATURE=%g\n'',p(4));', ...
[ 'fprintf(fid, ''# number of ion types and masses\n  NTYPES = ' num2str(numel(geom.symbols)) '\n  MASS = ' mass '\n'');' ], ...  
  'fprintf(fid, ''# q points section\n  LRECIP = .T.\n'');', ...
  'sz0 = size(t);', ...
  'if all(cellfun(@(c)max(size(c)) == numel(c), {x y z})) && numel(unique(cellfun(@(x)length(x), {x y z}))) > 1', ...
  '  fprintf(fid, ''  ND=%i\n'', numel(x));', ...
  '  fprintf(fid, ''  NPOINTS=1\n'');', ...
  '  % first write the min locations QI=QF', ...
  '  fprintf(fid, ''  QI ='');', ...
  '  for i1=1:numel(x)', ...
  '    fprintf(fid, ''    %f %f %f'', x(i1),y(i1),z(i1));', ...
  '    if i1<numel(x), fprintf(fid,'' \\\n'');', ...
  '    else fprintf(fid,''\n''); end', ...
  '  end % for QI', ...
  '  % then write the max locations QF=QI', ...
  '  fprintf(fid, ''  QF ='');', ...
  '  for i1=1:numel(x)', ...
  '    fprintf(fid, ''    %f %f %f'', x(i1),y(i1),z(i1));', ...
  '    if i1<numel(x), fprintf(fid,'' \\\n'');', ...
  '    else fprintf(fid,''\n''); end', ...
  '  end % for QF', ...
  'else', ...
  '  xu=unique(x(:)); yu=unique(y(:)); zu=unique(z(:));', ...
  '  sz = [ numel(xu) numel(yu) numel(zu) ];', ...
  '  [dummy,L_id] = max(sz);             % the dimension which will be a line ND', ...
  '  G_id = find(L_id ~= 1:length(sz));  % the other 2 smaller sizes in grid', ...
  '  fprintf(fid, ''  ND=%i\n'', prod(sz(G_id)));', ...
  '  fprintf(fid, ''  NPOINTS = %i\n'', sz(L_id));', ...
  '  % first write the min locations', ...
  '  fprintf(fid, ''  QI ='');', ...
  '  for i1=1:sz(G_id(1))', ...
  '    for i2=1:sz(G_id(2))', ...
  '      if L_id == 1, X=min(xu); Y=yu(i1); Z=zu(i2);', ...
  '      elseif L_id == 2, Y=min(yu); X=xu(i1); Z=zu(i2);', ...
  '      else Z=min(zu); X=xu(i1); Y=yu(i2); end', ...
  '      fprintf(fid, ''    %f %f %f'', X,Y,Z);', ...
  '      if i1*i2~=prod(sz(G_id)), fprintf(fid,'' \\\n'');', ...
  '      else fprintf(fid,''\n''); end', ...
  '    end', ...
  '  end % for QI', ...
  '  % then write the max locations', ...
  '  fprintf(fid, ''  QF ='');', ...
  '  for i1=1:sz(G_id(1))', ...
  '    for i2=1:sz(G_id(2))', ...
  '      if L_id == 1, X=max(xu); Y=yu(i1); Z=zu(i2);', ...
  '      elseif L_id == 2, Y=max(yu); X=xu(i1); Z=zu(i2);', ...
  '      else Z=max(zu); X=xu(i1); Y=yu(i2); end', ...
  '      fprintf(fid, ''    %f %f %f'', X,Y,Z);', ...
  '      if i1*i2~=prod(sz(G_id)), fprintf(fid,'' \\\n'');', ...
  '      else fprintf(fid,''\n''); end', ...
  '    end', ...
  '  end % for QF', ...
  'end % 3D grid', ...
  'fclose(fid);', ...
  'clear xu yu zu X Y Z', ...
  '% call PHON for FREQ', ...
  'try', ...
  '  cd(target);', ...
[ '  [status,result] = system(''' precmd options.phon ' > phon.log'');' ], ...
  '  % import FREQ', ...
  '  FREQ=load(''FREQ'',''-ascii'')/.24180; % THz -> meV', ...
  'catch; disp(fileread(''phon.log'')); ', ...
  '  disp([ ''model '' this.Name '' '' this.Tag '' could not run PHON from '' target ]); return', ...
  'end', ...
  'save FREQ.mat FREQ', ...
  'try; delete(''FREQ''); end', ...
  'try; delete(''FREQ.cm''); end', ...
  'try; delete(''INPHON''); end', ...
  '% get the DOS if the grid is 3D', ...
  'if all(cellfun(@ndims,{x y z}) == 3)', ...
  '  DOS = load(''DOS'',''-ascii''); DOS=iData(DOS(:,1)/.24180,DOS(:,2)/sum(DOS(:,2)));', ...
  'else DOS=[];', ...
  'end', ...
  'if all(cellfun(@numel,{x y z}) > 9)  % when not quick HKLE computation', ...
  '  % compute DOS on a regular grid and thermodynamic properties', ...
  '  % write the INPHON for DOS and THERMO -------------------------------------', ...
[ '  fid=fopen(fullfile(target,''INPHON''),''w'');' ], ...
[ '  fprintf(fid, ''# Control file for PHON to compute the TERMO/DOS from POSCAR and FORCES in %s\n'', target);' ], ...
  '  fprintf(fid, ''#   PHON: <http://chianti.geol.ucl.ac.uk/~dario>\n'');', ...
  '  fprintf(fid, ''#   D. Alfe, Computer Physics Communications 180,2622-2633 (2009)\n'');', ...
  '  fprintf(fid, ''  LSYMM=.TRUE.\n'');', ...
  '  fprintf(fid, ''  LSUPER=.F.\n'');', ...
  '  fprintf(fid, ''  LFREE=.T.\n'');', ...
  '  fprintf(fid, ''  TEMPERATURE=%g\n'',p(4));', ...
[ '  fprintf(fid, ''# NDIM   = ' num2str(options.supercell) '\n'');' ], ...
[ '  fprintf(fid, ''# number of ion types and masses\n  NTYPES = ' num2str(numel(geom.symbols)) '\n  MASS = ' mass '\n'');' ], ...  
  '  fprintf(fid, ''  DOSSMEAR = 0.02\n'');', ...
  '  fprintf(fid, ''  QA=11; QB=11; QC=11\n'');', ...
  '  fclose(fid);', ...
  '  try', ...
[ '    [status,result] = system(''' precmd options.phon ' >> phon.log'');' ], ...
  '    if isempty(DOS)', ...
  '      DOS = load(''DOS'',''-ascii'');', ...
  '      if p(5) > 0, DOS(:,1)=DOS(:,1)*p(5); end', ...
  '      DOS=iData(DOS(:,1)/.24180,DOS(:,2)/sum(DOS(:,2)));', ...
  '    end', ...
  '    % get THERMO properties', ...
  '    fid=fopen(''THERMO'');', ...
  '    if fid ~= -1', ...
  '      THERMO=textscan(fid,''%f'',''CommentStyle'',''#''); THERMO=THERMO{1}; fclose(fid);', ...
  '      this.UserData.properties.temperature=THERMO(1);', ...
  '      this.UserData.properties.zero_point_energy=THERMO(2);', ...
  '      this.UserData.properties.free_energy=THERMO(3);', ...
  '      this.UserData.properties.internal_energy=THERMO(4);', ...
  '      this.UserData.properties.entropy=THERMO(5);', ...
  '      this.UserData.properties.specific_heat=THERMO(6);', ...
  '    end', ...
  '  catch; ', ...
  '    disp([ ''model '' this.Name '' '' this.Tag '' could not get THERMO PHON from '' target ]);', ...
  '  end', ...
  'end', ...
  'cd(pw);', ...
  'if ~isempty(DOS)', ...
  '  DOS.Title = [ ''DOS '' this.Name ]; xlabel(DOS,''DOS Energy [meV]''); DOS.Error=0;', ...
  '  this.UserData.DOS=DOS;', ...
  'end', ...
  '% multiply all frequencies(columns, meV) by a DHO/meV', ...
  'Amplitude = p(1); Gamma=p(2); Bkg = p(3); T=p(4);', ...
  'if p(5) > 0, FREQ=FREQ*p(5); end', ...
  script_dho{:} };

signal=iFunc(signal);

% ------------------------------------------------------------------------------

function [poscar, options]=sqw_phon_argin(poscar, options)

  options.calculator = 'QuantumEspresso';

  % random displacement
  if isfield(options,'disp') && strcmp(options.disp,'random')
    % the DISP initial vector is set with a random set of [-1 0 1]
    options.disp = round(2*rand(1,3)-1);
    options.disp(abs(options.disp) > 3) = sign(options.disp(abs(options.disp) > 3));
    this = find(~options.disp); if numel(this) == 2, options.disp=sign(randn); end
  end

  % copy the POSCAR into the target directory
  d = dir(poscar);
  try
    copyfile(poscar, options.target);
  end
  % copy any additional potential
  try
    copyfile(fullfile(fileparts(poscar),'*.UPF'),options.target);
  end

  poscar = fullfile(options.target, d.name);

% ------------------------------------------------------------------------------

function options = sqw_phon_potentials(geom, options)

% missing potentials: use QE default location
if numel(options.potentials) < numel(geom.atomcount)
  options.potentials{end+1} = fullfile(getenv('HOME'),'espresso','pseudo');
  options.potentials{end+1} = getenv('PSEUDO_DIR');
  options.potentials{end+1} = getenv('ESPRESSO_PSEUDO');
  options.potentials{end+1} = pwd;
  options.potentials{end+1} = fileparts(geom.filename);
  if isunix
    options.potentials{end+1} = '/usr/share/espresso/pseudo';
    options.potentials{end+1} = '/usr/local/espresso/pseudo';
  end
end

% potential elements given as path: we get all potentials there in
potentials = {}; potentials_full = {};
if ischar(options.potentials), options.potentials = { options.potentials }; end
for index=1:numel(options.potentials)
  if isdir(options.potentials{index})
    % add the full content of the directory
    d = dir(options.potentials{index});
    for i=1:numel(d)
      if ~d(i).isdir && ~any(strcmp(d(i).name, potentials))
        [p,f,e] = fileparts(fullfile(options.potentials{index},d(i).name));
        if strcmp(lower(e),'.upf')
          potentials_full{end+1} = fullfile(options.potentials{index},d(i).name);
          potentials{end+1}      = d(i).name; 
        end
      end
    end
  else
    this = dir(options.potentials{index});
    if ~isempty(this) && ~any(strcmp(this.name, potentials)) % must exist
      potentials_full{end+1} = options.potentials{index};
      potentials{end+1}      = this.name; 
    end
  end
end
options.potentials      = potentials;
options.potentials_full = potentials_full;

% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
function geom1 = sqw_phon_supercell(poscar, options)
% sqw_phon_supercell: create a INPHON file for PHON to create a supercell POSCAR
%   returns the number of atoms, number of atom types, name... in a structure
%
% http://www.mathworks.com/matlabcentral/fileexchange/36836-vasplab/content//vasplab/import_poscar.m then supercell then export_poscar

if isempty(options.supercell),  options.supercell=2; end
if isscalar(options.supercell), options.supercell=[ options.supercell options.supercell options.supercell ]; end

geom1 = import_poscar(poscar);

% analyse chemical formula
if isempty(geom1.symbols)
  sqw_phon_error([ mfilename ': the file ' poscar ' MUST contain the chemical formula of the compound (1st line or before number of atoms line).' ], options);
  geom1=[]; return
else
  % build the structure holding all elements with counts
  geom1.elements = [];
  symbols = {};
  for index=1:numel(geom1.symbols)
    try
      this   = parse_formula(geom1.symbols{index});
    catch
      continue
    end
    fields = fieldnames(this);
    for i=1:numel(fieldnames(this))
      if isfield(geom1.elements, fields{i})
        geom1.elements.(fields{i}) = geom1.elements.(fields{i}) + this.(fields{i});
      else
        geom1.elements.(fields{i}) = this.(fields{i});
      end
      if ~any(strcmp(fields{i}, symbols))
        symbols{end+1} = fields{i};
      end
    end
    
  end
  geom1.symbols = symbols;
  % the number of fields in elements should match the number of atom types
  if numel(fieldnames(geom1.elements)) ~= numel(geom1.atomcount)
    disp([ mfilename ': WARNING: the number of elements in the formula (' ...
      num2str(numel(fieldnames(geom1.elements))) ') does not match the number of '
      'atom counts (' numel(geom1.atomcount) ')' ]);
    geom1.elements
    geom1.atomcount
  end
end

% if the FORCES file was given at start, we use the POSCAR as a supercell
if isfield(options,'force_file')
  return
end

[p,f] = fileparts(poscar);
if isempty(p), p=pwd; end
% check if this is already a super cell or FORCES are here
if isempty(strfind(lower(geom1.comment),'supercell'))
  disp([ mfilename ': generating supercell and initial displacement guess from ' poscar ]);
  % generate INPHON for supercell generation
  fid = fopen(fullfile(p,'INPHON'),'w');
  if fid < 0
    sqw_phon_error([ mfilename ': could not create file ' fullfile(p,'INPHON') ], options);
    geom1=[]; return
  end
  fprintf(fid, '# Control file for PHON to generate a supercell from %s\n', poscar);
  fprintf(fid, '#   PHON: <http://chianti.geol.ucl.ac.uk/~dario>\n');
  fprintf(fid, '#   D. Alfe, Computer Physics Communications 180,2622-2633 (2009)\n');
  fprintf(fid, 'LSUPER =.TRUE.\n');
  fprintf(fid, 'NDIM   = %i %i %i\n', options.supercell);
  fprintf(fid, 'NTYPES = %i\n', numel(geom1.atomcount));
  fprintf(fid, 'IPRINT = 2\n');

  if isfield(options,'accuracy') && strcmp(options.accuracy,'fast')
    fprintf(fid, 'LCENTRAL = .F.\n');
  end

  if isfield(options,'disp')    
    if numel(options.disp) == 3 % given as a direction vector
      fprintf(fid, 'DXSTART= %f %f %f\n', int16(options.disp));
      disp([ mfilename ': using initial displacement [ ' num2str(int16(options.disp(:)')) ' ]' ]);
      options.disp = 0.01*norm(options.disp);
    end
    if isscalar(options.disp)
      % convert from Angs to PHON value: 25 -> 0.04 Angs
      options.disp = max(1/options.disp, 1);

      fprintf(fid, 'DISP   = %i\n', round(options.disp));
    end
  end
  fclose(fid);
  
  % make sure we use a POSCAR file
  if ~strcmp(f,'POSCAR')
    copyfile(poscar, fullfile(p,'POSCAR'),'f');
  end
  % required to avoid Matlab to use its own libraries
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ; DISPLAY= ; ';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; DISPLAY= ; '; 
  else           precmd=''; end
  
  % call phon in path 'p'. It also provides DISP
  if ~isempty(options.phon)
    pw = pwd;
    cd(p);
    disp([ 'cd(''' p '''); ' options.phon ]);
    [status, result] = system([ precmd options.phon ' > phon.log' ]);
    cd(pw);
  end
  % check if the expected files have been created
  if isempty(dir(fullfile(p,'SPOSCAR'))) || isempty(dir(fullfile(p,'DISP')))
    try; disp(fileread(fullfile(p,'phon.log'))); end
    sqw_phon_error([ mfilename ': Error executing PHON: could not create SPOSCAR and DISP file.' ], options);
    geom1=[]; return
  end
  % modify 1st line so that the initial system name is retained
  geom2 = import_poscar(fullfile(p,'SPOSCAR'));
  geom2.comment = [ strtrim(geom1.comment) ' supercell ' mat2str(options.supercell) ];
  geom2.symbols = geom1.symbols;
  geom2.elements= geom1.elements;
  
  if strcmp(f,'POSCAR') % make a copy before over-writing
    copyfile(poscar, [ poscar '_' datestr(now,30) ]);
  end
  % copy supercell into POSCAR
  export_poscar(fullfile(p,'POSCAR'), geom2);
  geom1 = geom2;
  geom1.filename = poscar;
else
  disp([ mfilename ': re-using supercell ' fullfile(p,'POSCAR') ]);
end
disp([ mfilename ': system: ' geom1.comment ]);

% ------------------------------------------------------------------------------
function [forces, geom] = sqw_phon_forces(geom, options)
% sqw_phon_forces: generate forces using QE

forces = [];
% check if FORCES are here
[p,f,e] = fileparts(geom.filename);
if isempty(p), p=pwd; end
if isempty(dir(fullfile(p,'FORCES')))
  disp([ mfilename ': generating Hellmann-Feynman FORCES from ' f e ]);
  % we read the DISP file
  % rows: [ n_at delta_x delta_y delta_z ]
  try
    displacements = fileread(fullfile(p,'DISP'));
  catch
    sqw_phon_error([ mfilename ': the displacement file ' fullfile(p,'DISP') ' is missing. Reset initial POSCAR file and re-run (not a supercell).' ], options);
    return
  end
  displacements(displacements == '"' | displacements == '\') = '';
  displacements = str2num(displacements);
  
  % find suitable potentials
  if ~isempty(options.command)
    [geom.potentials,geom.potentials_full] = sqw_phon_forces_pwscf_potentials(geom, options);
  end
  
  % determine the type of the atoms
  geom.type = [];
  for j=1:numel(geom.atomcount)
    geom.type = [ geom.type ones(1,geom.atomcount(j)) ];
  end
  
  options.status = sprintf('Starting computation. %i displacements.\n', size(displacements,1));
  sqw_phonons_htmlreport('', 'status', options);
  
  t0 = clock;           % a vector used to compute elapsed/remaining seconds
  
  % now for each displacement line, we move the atom in geom.coords
  % and then call QE or VASP
  
  for move=1:size(displacements,1)
    displaced = geom; % copy initial geometry (POSCAR)
    index = displacements(move, 1); % this is the atom ID which is to be moved in [1:natoms]
    if (index< 0 || index > size(geom.coords,1))
      disp([ mfilename ': displacement ' num2str(move) ' moves atom ' num2str(index) ]);
      disp(['  but there are only ' num2str(size(geom.coords,1)) ' coordinates in the supercell POSCAR.']);
      disp('   Skipping');
      continue;
    end
    
    % move the atom
    displaced.coords(index,:) = displaced.coords(index,:) ...
                              + displacements(move, 2:4); % move XYZ
    
    if ishandle(options.gui), waitbar(0.15+((move-1)/size(displacements,1)*.6), options.gui, [ mfilename ': moving atom ' num2str(move) '/' num2str(size(displacements,1)) ]); end
    disp([ mfilename ': step ' num2str(move) '/' num2str(size(displacements,1)) ...
        ' moving atom ' displaced.symbols{displaced.type(displacements(move, 1))} ...
        ' by ' mat2str(displacements(move, 2:4)) ]);
    
    % compute FORCES, then add the result to the 'forces' cell
    if ~isempty(options.command)
      forces{move}= sqw_phon_forces_pwscf(displaced, options);
    else
      forces{move} = [];
    end
    if isempty(forces{move})
      disp([ mfilename ': aborting FORCES computation in ' p ]);
      forces = [];
      return
    elseif move < size(displacements,1)
      % display ETA. There are size(displacements,1) steps.
      % up to now we have done 'move' and it took etime(clock, t)
      % time per iteration is etime(clock, t)/move.
      % total time of computation is etime(clock, t)/move*size(displacements,1)
      % time remaining is etime(clock, t)/move*(size(displacements,1)-move)
      % final time is     t+etime(clock, t)/move*size(displacements,1)
      remaining = etime(clock, t0)/move*(size(displacements,1)-move);
      hours     = floor(remaining/3600);
      minutes   = floor((remaining-hours*3600)/60);
      seconds   = floor(remaining-hours*3600-minutes*60);
      enddate   = addtodate(now, ceil(remaining), 'second');
      
      options.status = [ 'ETA ' sprintf('%i:%02i:%02i', hours, minutes, seconds) ', ending on ' datestr(enddate) ' [' num2str(round(move*100.0/size(displacements,1))) '%]' ];
      disp([ mfilename ': ' options.status ]);
      sqw_phonons_htmlreport('', 'status', options);
    end

  end %for index

  % WRITE the FORCES file
  fid = fopen(fullfile(p,'FORCES'), 'w');
  if fid < 0
    sqw_phon_error([ mfilename ': could not create file ' fullfile(p,'FORCES') ], options);
    return
  end
  fprintf(fid, '%i\n', numel(forces));  % nb of displacements
  for index=1:numel(forces)
    fprintf(fid, '%s\n', num2str(displacements(index,:)));
    this = cellstr(num2str(forces{index}));
    fprintf(fid, '    %s\n', this{:});
  end
  
else
  disp([ mfilename ': re-using ' fullfile(p,'FORCES') ]);
  disp('  Delete this file if you want to re-generate Hellman-Feynman FORCES e.g. with other pseudo-potentials.');
end

% ------------------------------------------------------------------------------
function force = sqw_phon_forces_pwscf(displaced, options)

  force = [];
  % check if FORCES are here
  [p,f] = fileparts(displaced.filename);
  if isempty(p), p=pwd; end

  % CONTROL: we write the control file for Quantum ESPRESSO 'pw.d'
  fid = fopen(fullfile(p,'pw.d'), 'w');
  if fid < 0
    sqw_phon_error([ mfilename ': could not create file ' fullfile(p,'pw.d') ], options);
    return
  end
  natoms = size(displaced.coords,1);
  ntyp   = numel(displaced.atomcount);
  alat   = norm(displaced.lattice(1,:))/0.529; % in Bohr [a.u] unit
  if isfield(options,'kpoints') && ~isempty(options.kpoints) && all(options.kpoints > 0)
    kpoints = options.kpoints; 
  else kpoints   = 2; 
  end
  if isscalar(kpoints), kpoints=[ kpoints kpoints kpoints ]; end
  if isfield(options,'ecutwfc') && ~isempty(options.ecutwfc) && options.ecutwfc > 0
    ecut = options.ecutwfc; 
  else ecut   = ntyp*15; 
  end
  if isfield(options,'mixing_beta') && ~isempty(options.mixing_beta) && options.mixing_beta > 0
    mixing_beta = options.mixing_beta; 
  else mixing_beta = 0.7; 
  end
  
  fprintf(fid, '%s\n', displaced.comment); % chemical formula/system
  fprintf(fid, 'crystal\n');
  fprintf(fid, '&control\n');
  fprintf(fid, '  calculation=''scf''\n');
  fprintf(fid, '  restart_mode=''from_scratch'',\n');
  fprintf(fid, '  prefix=''%s''\n', strtok(displaced.comment));
  fprintf(fid, '  pseudo_dir=''./''\n');
  fprintf(fid, '  tprnfor=.true.\n');
  fprintf(fid, '/\n');
  fprintf(fid, '&system\n');
  fprintf(fid, '  ibrav = 0, celldm(1) = %f\n', alat);
  fprintf(fid, '  nat= %i, ntyp= %i,\n', natoms, ntyp);
  fprintf(fid, '  ecutwfc = %f\n', ecut);
  if isfield(options,'ecutrho')
    fprintf(fid, '  ecutrho = %f\n', options.ecutrho);
  end
  if isfield(options,'nbnd')
    fprintf(fid, '  nbnd = %i\n', options.nbnd);
  end
  if isfield(options,'occupations') && ~isempty(options.occupations)
    if ischar(options.occupations)
      switch lower(options.occupations)
      case {'smearing','metal'}
      fprintf(fid, '  occupations=''smearing'', smearing=''methfessel-paxton'', degauss=0.04\n');
      case {'fixed','insulator'}
      fprintf(fid, '  occupations=''fixed''\n');
      otherwise
      fprintf(fid, '  occupations=''%s''\n', options.occupations);
      end
    elseif isscalar(options.occupations)
      fprintf(fid, '  occupations=''smearing'', smearing=''methfessel-paxton'', degauss=%g\n', ...
        options.occupations);
    end
  end
  fprintf(fid, '/\n');
  fprintf(fid, '&electrons\n');
  if isfield(options, 'toldfe') && ~isempty(options.toldfe)
    fprintf(fid, '    conv_thr = %f\n', options.toldfe);
  end
  fprintf(fid, '    mixing_beta = %f\n', mixing_beta);
  fprintf(fid, '    mixing_mode = ''plain''\n');
  if isfield(options,'mixing_ndim') && ~isempty(options.mixing_ndim)
    fprintf(fid, '    mixing_ndim = %i\n', options.mixing_ndim);
  end
  if isfield(options, 'electron_maxstep') && ~isempty(options.electron_maxstep)
    fprintf(fid, '    electron_maxstep = %i\n', options.electron_maxstep);
  end
  if isfield(options, 'diagonalization') && ~isempty(options.diagonalization)
    if strncmp(options.diagonalization, 'dav',3)
      options.diagonalization = 'david';
    elseif ~strncmp(options.diagonalization,'cg',2)
      options.diagonalization = 'cg';
    end
    fprintf(fid, '    diagonalization = ''%s''\n', options.diagonalization);
  end
  fprintf(fid, '/\n');
  fprintf(fid, 'ATOMIC_SPECIES\n');
  for index=1:numel(displaced.symbols)
  % write: Element Mass pseudo-potential-file
    fprintf(fid,'%3s %5f %s\n', displaced.symbols{index}, ...
      molweight(displaced.symbols{index}), displaced.potentials{index});
    % the potentials must be copied locally to be visible by PWSCF
    try
    copyfile(displaced.potentials_full{index}, p);
    end
  end
  fprintf(fid, 'ATOMIC_POSITIONS { crystal }\n');
  for index=1:size(displaced.coords,1)
  % write: Element x y z
    fprintf(fid, '%3s %f %f %f\n', displaced.symbols{displaced.type(index)}, displaced.coords(index,:));
  end

  fprintf(fid, 'K_POINTS { automatic }\n');
  fprintf(fid, ' %i %i %i 0 0 0\n', kpoints);
  fprintf(fid, 'CELL_PARAMETERS { alat }\n');
  % normalise the cell parameters
  for index=1:3, 
   displaced.lattice(index,:)=displaced.lattice(index,:)/norm(displaced.lattice(index,:)); 
  end
  s = cellstr(num2str(displaced.lattice)); % 
  fprintf(fid, '%s\n', s{:});
  fprintf(fid, '\n');
  fclose(fid);

  % EXEC: we run QE/pw.x and collect stdout
  % required to avoid Matlab to use its own libraries
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ; DISPLAY= ; ';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; DISPLAY= ; '; 
  else           precmd=''; end
  disp([ options.command ' < ' fullfile(p,'pw.d') ' > ' fullfile(p, 'pw.out') ]);
  pw = pwd;
  cd(p);

  if isfield(options, 'mpi') && ~isempty(options.mpi) && options.mpi > 1
    [status, result] = system([ precmd options.mpirun ' ' options.command ' < pw.d > pw.out' ]);
  else
    [status, result] = system([ precmd options.command ' < pw.d > pw.out' ]);
  end

  cd(pw);
  % clear wavefunctions, which can be very large
  delete(fullfile(p, '*.wfc*'));

  % READ the QE/PWSCF output file and search for string 'Forces acting'
  L = fileread(fullfile(p, 'pw.out'));
  
  ismetallic = strfind(L, 'the system is metallic');
  if ~isempty(ismetallic)
    sqw_phon_error([ mfilename ': The system is metallic but you have used occupations=''fixed''. Rebuild model with occupations=''smearing''.' ], options)
    disp('occupations=')
    disp(options.occupations);
    return
  end

  forces_acting = strfind(L, 'Forces acting');
  if isempty(forces_acting)
    disp(L);
    disp([ mfilename ': convergence NOT achieved.' ]);
    disp([ 'TRY: sqw_phonons(..., ''mixing_beta=0.3; nsteps=200; toldfe=1e-6; occupations=smearing; ecut=' num2str(round(ecut*1.5*13.6)) ''')' ])
    sqw_phon_error([ mfilename ': PWSCF convergence NOT achieved.' ], options)
    return
  end
  L = L(forces_acting:end);
  % then read next natoms lines starting with 'atom', e.g.:
  % Forces acting on atoms (Ry/au):
  %
  % atom    1 type  1   force =    -0.00456135   -0.00456135   -0.00000000
  % ...
  atom = strfind(L, 'atom ');
  for index=1:numel(atom)
    this = L(atom(index):(atom(index)+80));
    V = sscanf(strtrim(this), 'atom %i type %i force = %f %f %f');

    if numel(V) == 5
      V = V(:)';
      if isempty(force), force = V(3:5); 
      else force = [ force ; V(3:5) ]; end
    end
    if size(force, 1) >= natoms, break; end % exit when we have one per atom in supercell
  end

  force = force *25.711; % from Ry/a.u to eV/Angs
  
% ------------------------------------------------------------------------------
  
function [potentials, potentials_full] = sqw_phon_forces_pwscf_potentials(displaced, options)

  % search for a suitable potential
  potentials = {}; potentials_full = {};
  for index=1:numel(displaced.symbols)
    % identify element in the list of potentials
    match = strcmpi(displaced.symbols{index}, strtok(options.potentials,'._'));
    match_full = options.potentials_full(match); match = options.potentials(match);
    if isempty(match)
      disp([ mfilename ': no suitable pseudo-potential found for atom ' displaced.symbols{index} ])
      disp('  Get such potentials at <http://www.quantum-espresso.org/pseudopotentials/>')
      disp([ '  and copy them in ' pwd ' or ' fullfile(getenv('HOME'),'espresso','pseudo') ]);
      if isempty(getenv('PSEUDO_DIR'))
        disp('  or define environment variable PSEUDO_DIR, and put PP in this location.');
      else
        disp([ '  or in ' getenv('PSEUDO_DIR') ]);
      end
      if ~isempty(options.potentials)
        disp('Available potentials are:')
        disp(options.potentials_full(:));
      end
      sqw_phon_error([ mfilename ': pseudo-potential missing for atom ' displaced.symbols{index} ], options);
      return
    elseif numel(match) > 1
      % more than one match: select PBE-PAW if possible
      pbe = find(~cellfun('isempty', strfind(lower(match), 'pbe')));
      paw = find(~cellfun('isempty', strfind(lower(match), 'paw')));
      if ~isempty(pbe)
        select = pbe(1);
      elseif ~isempty(paw)
        select = paw(1);
      else select=1;
      end
      % if still more than one choice, pop-up list selector
      if numel(match) > 1 && ishandle(options.gui)
        [select,OK] = listdlg('ListString', match, ...
            'ListSize', [400 200], ...
            'Name', [ 'Pseudo-potential for ' displaced.symbols{index} ], ...
            'InitialValue', select, ...
            'PromptString', { [ mfilename ': Select the pseudo-potential ', ...
                              'to use for atom ' displaced.symbols{index} ], ...
                              'The PBE-PAW potentials are recommended.', ...
                              'Get more pseudo-potentials at:', ...
                              'http://www.quantum-espresso.org/pseudopotentials/' });
        if isempty(select), return; end
        match = match(select); match_full = match_full(select);
      end

    end
    potentials{index} = match{1}; potentials_full{index} = match_full{1};
  end
  

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

function sqw_phon_error(message, options)

if ishandle(options.gui)
  delete(options.gui);
  errordlg(message, [ 'iFit: ' mfilename ' ' options.configuration ' FAILED' ]);
end
if ~isdeployed && usejava('jvm') && usejava('desktop')
  disp([ '<a href="matlab:doc(''sqw_phonons'')">sqw_phonons help</a>' ])
end
disp(message);
