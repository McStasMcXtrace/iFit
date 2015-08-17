function signal=sqw_phon(varargin)
% model=sqw_phon(poscar, pseudo, ..., options)  
%
%   iFunc/sqw_phon: computes phonon dispersions from Ab-Initio.
%   A model which compute phonon dispersions from the forces acting between
%     atoms. The input arguments is a VASP-type POSCAR file providing the
%     geometry of the cell (lattice and atom positions). 
%     It is recommended to specify if the system is a metal or an insulator.
%     Additional arguments can be given to control the procedure used to
%     compute a supercell and the forces using QuantumEspresso (ab-initio).
%
% WARNING: Single intensity and line width parameters are used here.
%   This model is only suitable to compute phonon dispersions for e.g solid-
%   state materials.
%
% The arguments can be any set of:
% POSCAR: file name to an existing POSCAR
%   a VASP compliant file providing initial system structure, with lattice
%   and position of atoms in a stable, equilibrated configuration.
%   The first line MUST start with the chemical formula of the system, as a 
%   single word, which can be followed by any other comment.
%   The specified file name can be different than POSCAR.
%
% pseudo: file name of pseudo-potential to use. If not set, a reasonable 
%   choice will be made from existing UPF files from Quantum Espresso, with
%   a preference for PBE-PAW pseudo-potentials. Potentials are searched in the
%   QE location, in PSEUDO_DIR environment variable (when set) and locally. 
%   You can specify as many of these UPF and directories.
%
% 'metal' or 'insulator': indicates the type of occupation for electronic states.
%
% 'random': indicates that the initial displacement should be chosen randomly.
%
% options: a structure with optional settings:
%   options.target =path                   where to store all files and FORCES
%     a temporary directory is created when not specified.
%   options.NDIM   =scalar or [nx ny nz]   supercell size
%   options.disp   =[dx dy dz]             initial displacement direction
%     the initial displacement can e.g. be [ 1 -1 1 ]
%   options.potentials={ cell of strings } list of UPF potentials to use
%     if not given a PBE-PAW potential will be used (when available)
%     these can also be directories, which content is then added to the list.
%   options.kpoints=scalar or [nx ny nz]   Monkhorst-Pack grid
%   options.mpi    =scalar                 number of CPUs to use for PWSCF
%     this option requires MPI to be installed (e.g. openmpi).
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
%     wavefunctions. default is 15*ntyp in [Ry]. Larger value improves convergence.
%   options.electron_maxstep=scalar        max number of iterations for SCF.
%     default=100. Larger value improves convergence.
%   options.conv_thr=scalar                Convergence threshold for 
%     selfconsistency. default=1e-6.
%
% The options can also be entered as a single string with 'field=value; ...'.
%
% Once the model has been created, its use requires that axes are given on
% regular qx,qy,qz grids.
%     
% Example:
%   s=sqw_phon([ ifitpath 'Data/POSCAR_Al'],'metal','mpi=4');
%   qh=linspace(0,.5,50);qk=qh; ql=qh; w=linspace(0.01,100,51);
%   f=iData(s,[],qh,qk,ql,w); scatter3(log(f(1,:, :,:)),'filled');
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
%          or p='guess'
%         qh: axis along QH in rlu (row,double)
%         qk: axis along QK in rlu (column,double)
%         ql: axis along QL in rlu (page,double)
%         w:  axis along energy in meV (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value

% TODO: 
%   Write POSCAR/FORCES in tempname at eval from UserData if not available.

signal = [];

% ********* handle input arguments *********
[poscar, options]=sqw_phon_argin(varargin{:});

% check for installation of PHON and QE
[options.phon, options.pwscf] = sqw_phon_requirements;
if isempty(options.phon) || isempty(options.pwscf)
  return
end

% ********* compute the forces *********

% check/create supercell
geom = sqw_phon_supercell(poscar, options); % 2x2x2 supercell

% ********* look for potentials *********

options = sqw_phon_potentials(geom, options);

% check for FORCES
[forces, geom] = sqw_phon_forces(geom, options);

mass = [];
compound = ''; potentials = '';
for index=1:numel(geom.symbols)
  mass(index) = molweight(geom.symbols{index});
  compound    = [ compound geom.symbols{index} num2str(geom.atomcount(index)) ' ' ];
  potentials  = [ potentials geom.potentials{index} ' ' ];
end
mass = num2str(mass);

% ************** BUILD THE MODEL **************
signal.Name           = [ compound 'S(q,w) 3D dispersion PHON/QuantumEspresso with DHO line shape [' mfilename ']' ];

signal.Description    = [ compound 'S(q,w) 3D dispersion PHON/QuantumEspresso with DHO line shape. Pseudo-Potentials: ' potentials ];

signal.Parameters     = {  ...
  'Amplitude' ...
  'Gamma Damped Harmonic Oscillator width in energy [meV]' ...
  'Background' ...
  'Temperature [K]' ...
   };
  
signal.Dimension      = 4;         % dimensionality of input space (axes) and result

signal.Guess = [ 1 .1 0 10 ];

p = fileparts(poscar);

signal.UserData.POSCAR = fileread(poscar);
signal.UserData.FORCES = fileread(fullfile(p,'FORCES'));
signal.UserData.dir    = p;
signal.UserData.options= options;
signal.UserData.potentials=potentials;

signal.Expression     = { ...
  '% check if FORCES and POSCAR are here', ...
[ 'if isempty(dir(''' poscar ''')) || isempty(dir(''' fullfile(p,'FORCES') '''))' ], ...
[ '  error([ ''model '' this.Name '' '' this.Tag '' requires POSCAR and FORCES to reside in ' p ''' ])' ], ...
  'end', ...
  '% write INPHON for FREQ generation', ...
[ '  pw = pwd; cd(''' p ''');' ], ...
  '  % write the INPHON', ...
[ '  fid=fopen(''INPHON'',''w'');' ], ...
[ '  fprintf(fid, ''# Control file for PHON to compute the modes from POSCAR and FORCES in ' p '\n'');' ], ...
  '  fprintf(fid, ''#   PHON: <http://chianti.geol.ucl.ac.uk/~dario>\n'');', ...
  '  fprintf(fid, ''#   D. Alfe, Computer Physics Communications 180,2622-2633 (2009)\n'');', ...
  '  fprintf(fid, ''  LSYMM=.TRUE.\n'');', ...
  '  fprintf(fid, ''  LSUPER=.F.\n'');', ...
[ '  fprintf(fid, ''# number of ion types and masses\n  NTYPES = ' num2str(numel(geom.symbols)) '\n  MASS = ' mass '\n'');' ], ...  
  '  fprintf(fid, ''# q points section\n  LRECIP = .T.\n'');', ...
  '  sz0 = size(t);', ...
  '  if ndims(x) == 4, x=squeeze(x(:,:,:,1)); y=squeeze(y(:,:,:,1)); z=squeeze(z(:,:,:,1)); t=squeeze(t(1,1,1,:)); end',...
  'if max(size(x)) == numel(x) && max(size(y)) == numel(y) && max(size(z)) == numel(z)', ...
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
  '  fclose(fid);', ...
  '  % call PHON', ...
  'try', ...
[ '  [status,result] = system(''' options.phon ' > phon.log'');' ], ...
  '  % import FREQ', ...
  '  FREQ=load(''FREQ'',''-ascii'')/.24180; % THz -> meV', ...
  '  delete(''FREQ.cm'');', ...
  'end', ...
  '  cd(pw);', ...
  '  % multiply all frequencies(columns, meV) by a DHO/meV', ...
  '  Amplitude = p(1); Gamma=p(2); Bkg = p(3); T=p(4);', ...
  '  if T<=0, T=300; end', ...
  '  w=t(:) * ones(1,size(FREQ,1));', ...
  '  signal=zeros(size(w));', ...
  'for index=1:size(FREQ,2)', ...
  '% transform w and w0 to same size', ...
  '  w0= ones(numel(t),1) * FREQ(:,index)'';', ...
  '  toadd = Amplitude*Gamma *w0.^2.* (1+1./(exp(abs(w)/T)-1)) ./ ((w.^2-w0.^2).^2+(Gamma*w).^2);', ...
  '  signal = signal +toadd;', ...
  'end', ...
  'signal = reshape(signal'',sz0);'
};

signal=iFunc(signal);

%  '  if T<=0, T=300; end', ...
%  '  %p=ones(size(FREQ,2),1) * [ p(1) 0 p(2) p(3) p(4)/11.609 ];',...
%  '  %signal=zeros(size(p));', ...
%  '  %w=reshape(t,size(signal));', ...
%  '  %for index=1:size(FREQ,2)', ...
%  '  %  p(:,2) = FREQ(:,index);', ...
%  '  %  signal = signal+ p(3) *p(2)^2.*(1+1./(exp(abs(w)/p(4))-1))./((w.^2-p(2)^2).^2+(p(3)*t).^2);
%  '  %end', ...
%  '  %signal = signal*p(1);',...

% at model evaluation:


% call PHON

% read FREQ (THz)

% apply DHO to all energies

% when model is successfully built, display citations for PHON and QE
disp([ 'Model ' compound ' built using: (please cite)' ])
disp(' * PHON: D. Alfe, Computer Physics Communications 180,2622-2633 (2009)')
disp('     <http://chianti.geol.ucl.ac.uk/~dario>. BSD license.')
disp(' * Quantum Espresso: P. Giannozzi, et al J.Phys.:Condens.Matter, 21, 395502 (2009)')
disp('     <http://www.quantum-espresso.org/>. GPL2 license.')
disp(' * iFit: E. Farhi et al, J. Neut. Res., 17 (2013) 5.')
disp('     <http://ifit.mccode.org>');

% ------------------------------------------------------------------------------

function [poscar, options]=sqw_phon_argin(varargin)

  poscar = [];
  options.NDIM   =2;
  options.potentials ={};
  options.kpoints=2;
  options.disp   =[];

  % read input arguments
  for index=1:numel(varargin)
  if ischar(varargin{index}) && isempty(dir(varargin{index}))
    % first try to build a structure from the string
    this = str2struct(varargin{index});
    if isstruct(this)
      varargin{index} = this;
    end
  end
  if ischar(varargin{index})
    [p,f,e] = fileparts(varargin{index});
    % handle static options: metal,insulator, random
    if strcmp(varargin{index},'smearing') || strcmp(varargin{index},'metal')
      options.occupations = 'smearing';
    elseif strcmp(varargin{index},'fixed') || strcmp(varargin{index},'insulator')
      options.occupations = 'fixed';
    elseif strcmp(varargin{index},'random')
      options.disp='random';
    elseif ~isempty(dir(varargin{index})) && ~isdir(varargin{index}) ...
        && (isempty(e) || ~strcmp(lower(e),'.upf'))
      % found an existing file, not a pseudo potential
      if isempty(poscar)
        poscar = varargin{index};
      end
    elseif strcmp(lower(e),'.upf') || isdir(varargin{index})
      % found a '.upf' file or directory: pseudo-potential
      options.potentials{end+1} = varargin{index};
    end
  elseif isstruct(varargin{index})
    % a structure: we copy the fields into options.
    this = varargin{index};
    f    =fieldnames(this);
    for i=1:numel(fieldnames(this))
      options.(f{i}) = this.(f{i});
    end
  end
  
end

% random displacement
if strcmp(options.disp,'random')
  % the DISP initial vector is set with a random set of [-1 0 1]
  options.disp = round(randn(1,3));
  options.disp(abs(options.disp) > 3) = sign(options.disp(abs(options.disp) > 3));
  this = find(~options.disp); if numel(this) == 2, options.disp=sign(randn); end
end

if isempty(poscar)
  poscar = 'POSCAR';
end

if ~isfield(options,'target')
  options.target = tempname; % everything will go there
  mkdir(options.target)
end
;
% copy the POSCAR into the target directory
d = dir(poscar);
copyfile(poscar, options.target);
% copy any additional potential
try
copyfile(fullfile(fileparts(poscar),'*.UPF'),options.target);
end

poscar = fullfile(options.target, d.name);
disp([ mfilename ': copying initial ' d.name ' into ' options.target ]);


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

function [phon, pwscf]  = sqw_phon_requirements
% check if PHON and QuantumExpresso are installed

% check for PHON
cmd = 'phon';
if ispc, cmd=[ cmd '.exe' ]; end
if ~isempty(dir('INPHON')), delete('INPHON'); end
[status, result] = system(cmd);
if isempty(strfind(upper(result),'PHON, VERSION'))
  phon = [];
  disp([ mfilename ': ERROR: requires PHON to be installed.' ])
  disp('  1- Get it at <http://chianti.geol.ucl.ac.uk/~dario>.');
  disp('  2- Extract, go in directory, compile with: ./configure; make;')
  disp('  3- copy executable src/phon into e.g. /usr/local/bin or /usr/bin')
  status = 'error: PHON not installed';
else
  phon = cmd;
end

% check for Quantum ESPRESSO Plane Wave Self Consistent Field
if ispc
  cmd = 'pw.exe';
else
  cmd = 'pw.x';
end
[status, result] = system([ 'echo 0 | ' cmd ]);
try
delete('input_tmp.in')
end
if isempty(strfind(upper(result),'PWSCF'))
  pwscf = [];
  disp([ mfilename ': ERROR: requires Quantum ESPRESSO to be installed.' ])
  disp('  1- Get it at <http://www.quantum-espresso.org>.');
  disp('  2- Install with e.g.:')
  if ~ispc
  disp('       * sudo apt-get install quantum-espresso (Debian class systems)')
  disp('       * sudo yum install quantum-espresso (RedHat class systems)')
  disp('       * compile from .tar.gz (other Linux/unix, MacOSX): ./configure; make pw; make install');
  else
  disp('       * double click e.g. qe-5.2.0-64bit-serial.exe')
  end
  status = 'error: Quantum ESPRESSO not installed';
else
  pwscf = cmd;
  delete('CRASH'); % this is created when nothing is done in QE/PWSCF
end

% ------------------------------------------------------------------------------
function geom1 = sqw_phon_supercell(poscar, options)
% sqw_phon_supercell: create a INPHON file for PHON to create a supercell POSCAR
%   returns the number of atoms, number of atom types, name... in a structure
%
% http://www.mathworks.com/matlabcentral/fileexchange/36836-vasplab/content//vasplab/import_poscar.m then supercell then export_poscar

if isempty(options.NDIM),  options.NDIM=2; end
if isscalar(options.NDIM), options.NDIM=[ options.NDIM options.NDIM options.NDIM ]; end

geom1 = import_poscar(poscar);

% analyse chemical formula
if isempty(geom1.symbols)
  error([ mfilename ': the file ' poscar ' MUST contain the chemical formula of the compound (1st line or before number of atoms line).' ]);
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

[p,f] = fileparts(poscar);
if isempty(p), p=pwd; end
% check if this is already a super cell
if isempty(strfind(lower(geom1.comment),'supercell'))
  disp([ mfilename ': generating supercell and initial displacement guess from ' poscar ]);
  % generate INPHON for supercell generation
  fid = fopen(fullfile(p,'INPHON'),'w');
  if fid < 0
    error([ mfilename ': could not create file ' fullfile(p,'INPHON') ]);
  end
  fprintf(fid, '# Control file for PHON to generate a supercell from %s\n', poscar);
  fprintf(fid, '#   PHON: <http://chianti.geol.ucl.ac.uk/~dario>\n');
  fprintf(fid, '#   D. Alfe, Computer Physics Communications 180,2622-2633 (2009)\n');
  fprintf(fid, 'LSUPER =.TRUE.\n');
  fprintf(fid, 'NDIM   = %i %i %i\n', options.NDIM);
  fprintf(fid, 'NTYPES = %i\n', numel(geom1.atomcount));

  if isfield(options,'disp')
    if sum(abs(options.disp)) >= 1
      fprintf(fid, 'DISP   = %i\n', round(25.0*sum(abs(options.disp))));  
    else
      fprintf(fid, 'DISP   = %i\n', 100);  
    end
    
    if numel(options.disp) == 3
      fprintf(fid, 'DXSTART= %f %f %f\n', options.disp);
      disp([ mfilename ': using initial displacement [ ' num2str(options.disp(:)') ' ]' ]);
      options.disp = options.disp/norm(options.disp);
    end
  end
  fclose(fid);
  
  % make sure we use a POSCAR file
  if ~strcmp(f,'POSCAR')
    copyfile(poscar, fullfile(p,'POSCAR'));
  end
  
  % call phon in path 'p'. It also provides DISP
  if ~isempty(options.phon)
    pw = pwd;
    cd(p);
    disp([ 'cd(''' p '''); ' options.phon ]);
    [status, result] = system(options.phon);
    cd(pw);
  end
  % check if the expected files have been created
  if isempty(dir(fullfile(p,'SPOSCAR'))) || isempty(dir(fullfile(p,'DISP')))
    error([ mfilename ': Error executing PHON: could not create SPOSCAR and DISP file ' ]);
  end
  % modify 1st line so that the initial system name is retained
  geom2 = import_poscar(fullfile(p,'SPOSCAR'));
  geom2.comment = [ strtrim(geom1.comment) ' supercell ' mat2str(options.NDIM) ];
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
    error([ mfilename ': the displacement file ' fullfile(p,'DISP') ' is missing. Reset initial POSCAR file and re-run (not a supercell).' ]);
  end
  displacements(displacements == '"' | displacements == '\') = '';
  displacements = str2num(displacements);
  
  % find suitable potentials
  if ~isempty(options.pwscf)
    [geom.potentials,geom.potentials_full] = sqw_phon_forces_pwscf_potentials(geom, options);
  end
  
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
    
    displaced.coords(index,:) = displaced.coords(index,:) ...
                              + displacements(move, 2:4); % move XYZ
    % determine the type of the atoms
    displaced.type = [];
    atomcount_id=1;
    atomcount = cumsum(geom.atomcount);
    for index=1:size(geom.coords,1)
      if index > atomcount(atomcount_id), atomcount_id=atomcount_id+1; end
      displaced.type(index)=atomcount_id;
    end
    if ~displaced.type
      disp([ mfilename ': can not determine the type of atom moved. Skipping.' ]);
      continue;
    end
    
    disp([ mfilename ': step ' num2str(move) '/' num2str(size(displacements,1)) ...
        ' moving atom ' displaced.symbols{displaced.type(displacements(move, 1))} ...
        ' by ' num2str(displacements(move, 2:4)) ]);
    
    % now use either QE or VASP
    if ~isempty(options.pwscf)
      this.displacement = [ index displacements(move,:) ];
      this.forces       = sqw_phon_forces_pwscf(displaced, options);
    else
      this.forces = [];
    end
    if isempty(this.forces)
      disp([ mfilename ': aborting FORCES computation.' ]);
      return
    end
  
    % add the result to the 'forces' matrix
    forces{move}      = this;

  end %for index

  % WRITE the FORCES file
  fid = fopen(fullfile(p,'FORCES'), 'w');
  if fid < 0
    error([ mfilename ': could not create file ' fullfile(p,'FORCES') ]);
  end
  fprintf(fid, '%i\n', numel(forces));  % nb of displacements
  for index=1:numel(forces)
    fprintf(fid, '%s\n', num2str(displacements(index,:)));
    this = cellstr(num2str(forces{index}.forces));
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
    error([ mfilename ': could not create file ' fullfile(p,'pw.d') ]);
  end
  natoms = size(displaced.coords,1);
  ntyp   = numel(displaced.atomcount);
  alat   = norm(displaced.lattice(1,:))/0.529; % in Bohr [a.u] unit
  if isfield(options,'kpoints') && options.kpoints > 0
    kpoints = options.kpoints; else kpoints   = 2; 
  end
  if isscalar(kpoints), kpoints=[ kpoints kpoints kpoints ]; end
  if isfield(options,'ecutwfc') && options.ecutwfc > 0
    ecut = options.ecutwfc; else ecut   = ntyp*15; 
  end
  if isfield(options,'mixing_beta') && options.mixing_beta > 0
    mixing_beta = options.mixing_beta; else mixing_beta = 0.7; 
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
  if isfield(options,'occupations')
    switch lower(options.occupations)
    case {'smearing','metal'}
    fprintf(fid, '  occupations=''smearing'', smearing=''methfessel-paxton'', degauss=0.04\n');
    case {'fixed','insulator'}
    fprintf(fid, '  occupations=''fixed''\n');
    otherwise
    fprintf(fid, '  occupations=''%s''\n', options.occupations);
    end
  end
  fprintf(fid, '/\n');
  fprintf(fid, '&electrons\n');
  if isfield(options, 'conv_thr')
    fprintf(fid, '    conv_thr = %f\n', options.conv_thr);
  end
  fprintf(fid, '    mixing_beta = %f\n', mixing_beta);
  fprintf(fid, '    mixing_mode = ''plain''\n');
  if isfield(options,'mixing_ndim')
    fprintf(fid, '    mixing_ndim = %i\n', options.mixing_ndim);
  end
  if isfield(options, 'electron_maxstep')
    fprintf(fid, '    electron_maxstep = %i\n', options.electron_maxstep);
  end
  if isfield(options, 'diagonalization')
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
  disp([ options.pwscf ' < ' fullfile(p,'pw.d') ' > ' fullfile(p, 'pw.out') ]);
  pw = pwd;
  cd(p)
  try
    if isfield(options, 'mpi')
      [status, result] = system([ 'mpirun -np ' num2str(options.mpi) ' ' options.pwscf ' < pw.d > pw.out' ]);
    else
      [status, result] = system([ options.pwscf ' < pw.d > pw.out' ]);
    end
  end
  cd(pw);

  % READ the QE/PWSCF output file and search for string 'Forces acting'
  L = fileread(fullfile(p, 'pw.out'));
  
  ismetallic = strfind(L, 'the system is metallic');
  if ~isempty(ismetallic)
    error([ mfilename ': The system is metallic but you have used occupations=''fixed''. Rebuild model with occupations=''smearing''.' ])
  end

  forces_acting = strfind(L, 'Forces acting');
  if isempty(forces_acting)
    disp([ mfilename ': convergence NOT achieved.' ]);
    disp([ 'TRY: sqw_phon(..., ''miximg_beta=0.3; electron_maxstep=200; conv_thr=1e-6; occupations=smearing; ecutwfc=' num2str(round(ecut*1.5)) ''')' ])
    error([ mfilename ': PWSCF convergence NOT achieved.' ])
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
    match = strcmpi(displaced.symbols{index}, strtok(options.potentials,'.'));
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
      error([ mfilename ': pseudo-potential missing for atom ' displaced.symbols{index} ]);
    elseif numel(match) > 1
      % more than one match: select PBE-PAW if possible
      pbe = find(~cellfun('isempty', strfind(lower(match), 'pbe')));
      paw = find(~cellfun('isempty', strfind(lower(match), 'paw')));
      if ~isempty(pbe)
        select = pbe(1);
      elseif ~isempty(paw)
        select = paw(1)
      else select=1;
      end
      % if still more than one choice, pop-up list selector
      if numel(match) > 1
        [select,OK] = listdlg('ListString', match, ...
            'ListSize', [400 200], ...
            'Name', [ 'Pseudo-potential for ' displaced.symbols{index} ], ...
            'InitialValue', select, ...
            'PromptString', { [ mfilename ': Select the pseudo-potential', ...
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
