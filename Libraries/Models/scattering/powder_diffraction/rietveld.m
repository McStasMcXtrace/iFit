function y = rietveld(varargin) 
% usage:
%   fits(rietveld, measurement, pars, options, constraints)
% with   rieltveld(instr_par, sample_par, instrument, mcstas_options, ...) as iFunc
%   and p={instr_par, sample_par} is a structure/cell which holds parameters that vary (numerics)
%          members as chars are assumed to remain constant.
% 
% Parameters of the model can be entered as:
% * a structure with one field 'structure.<atom>' per atom in the cell, named from the atom, e.g.
%     p.structure.Al1=[x y z {Biso occ spin charge}] where {Biso occ spin charge} are optional
%     p.Al1          =[...] can also be used for a compact specification (see exemple below).
%   plus additional fields
%     p.cell=[a b c alpha beta gamma] lattice cell parameters (Angs and deg)
%     p.Spgr='space group' as Number, Hall or Hermman-Mauguin
%
% * a structure to configure the McStas simulation parameters, with fields from
%       ncount: number of neutron events per iteration, e.g. 1e5 (double)
%       dir:    directory where to store results (string)
%       mpi:    number of processors/cores to use with MPI on localhost (integer) 
%       seed:   random number seed to use for each iteration (double)
%       gravitation: 0 or 1 to set gravitation handling in neutron propagation (boolean)
%       compile: 0 or 1 to force re-compilation of the instrument (boolean)
%       particle: 'n' (default) or 'x' depending if you use McStas or McXtrace (string)
%       monitors:  cell string of monitor names, or empty for all (cellstr)
%     Only the last monitor in the selection is used to compute the refinement criteria.
%     Type 'help mcstas' for more information about these items.
%
%   In addition, the 'CFML_write' field allows to set the name of the intermediate
%     reflection file written by the CrysFML routines, and usable by the PowderN
%     and Isotropic_Sqw McStas components. The default intermediate file is
%     'reflections.laz'.
%
% * a vector array where Z is the atomic number
%     p=[ a b c alpha beta gamma (Z x y z Biso occ spin charge) (...) ... ]
%
% * a CIF/CFL/PCR/shellX file name from which structure is extracted.
%
% * a McStas instrument description file name (.instr extension) to specify
%     the instrument model to use. The default instrument is 'templateDIFF.instr'
%
% * a char with label/value pairs separated by ';' to build a structure.
%
% * All additional structure fields are sent to the McStas instrument model.
%
% WARNING: at least ONE instrument parameter must be given.
%   All instrument parameters given as scalar numerical values are refined, 
%   others are kept fixed.
%
% Exemple: refine a NaCaAlF powder structure with the templateDIFF McStas instrument
%    Sample.title = 'Na2Ca3Al2F14';
%    Sample.cell  = [10.242696  10.242696  10.242696  90.000  90.000  90.000];
%    Sample.Spgr  = 'I 21 3';
%    Sample.Ca1   = [0.46737  0.00000  0.25000  0.60046  0.50000   0.0   2.0];
%    Sample.Al1   = [0.24994  0.24994  0.24994  0.18277  0.33333   0.0   3.0];
%    Sample.Na1   = [0.08468  0.08468  0.08468  2.25442  0.33330   0.0   1.0];
%    Sample.F1    = [0.13782  0.30639  0.12023  0.56669  1.00000   0.0  -1.0];
%    Sample.F2    = [0.36209  0.36363  0.18726  1.42140  1.00000   0.0  -1.0];
%    Sample.F3    = [0.46123  0.46123  0.46123  0.88899  0.33333   0.0  -1.0];
% The resulting hklF2 file 'reflections.laz' will be used by the instrument Powder 
% parameter and the model is computed at constant wavelength lambda=2.36 (given as char)
%    f = rietveld(Sample, 'templateDIFF.instr',' Powder=reflections.laz; lambda=''2.36'' ');
% Then 'f' is an iFunc Rietveld model. Then perform the Rietveld refinement.
%    p = fits(measurement, f); % where measurement holds a measured powder diffractogram
%
% See also: mcstas, iFunc, iData, iFunc/fits

% To initiate the model, one enters parameters as a structure with members.
% During evaluation, parameters are used as a vector matching initiated fields.

% extract parameters to build the CFL

y = [];

if nargin == 0
  disp([ mfilename ': require input arguments to define Rietveld model. See "help rietveld".'])
  disp('  Syntax: rietveld(structure,char,numerical,vector...). ')
  return
end

variable_p     = [];  % will store variable parameters for both Rietveld and instrument
constant_p     = [];  % will store fixed instrument parameters to be used by model
instrument     = '';  % name of the instrument model to use
Rietveld       = [];
mcstas_options = [];

atoms={'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar',...
      'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se',...
      'Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',...
      'Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',...
      'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au',...
      'Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu',...
      'Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt',...
      'Ds','Rg','Cn','Uut','Uuq','Uup','Uuh','Uuo'};

% get input parameters for the model ===========================================
% check input arguments:
% cell      -> treat iteratively
% char      -> convert to struct if appropriate
% char      -> read file with cif2hkl(--no-output-files --verbose <file>)
% numeric   -> [a b c aa bb cc (Z x y z B occ spin charge) .. ], use fixed instrument parameters
% structure -> analyse fields and search for Spgr, cell, atoms, x,y,z,...

for index=1:length(varargin)
  this = varargin{index};
  if ischar(this)   % char -> struct
    par = str2struct(this); 
    if isstruct(par), this = par; end
  end
  if ischar(this)                         % char -> file -> struct OR instrument
    % check extension
    [p,f,e] = fileparts(this);
    if strcmp(e, '.instr')
      instrument = this;
      continue; % jump to next argument in the list
    end
    % call cif2hkl
    if exist('cif2hkl') == 3
      % use MeX in verbose and no-output-files mode ('-')
      result = cif2hkl(this,[],[],'-',1);
      this = str2struct(result);
    else
      disp('cif2hkl is missing: compile it with e.g: ')
      disp('  gfortran -fPIC -c cif2hkl.f90')
      disp('  mex -O cif2hkl_mex.c cif2hkl.o -o cif2hkl -I/usr/lib/gcc/x86_64-linux-gnu/4.6 -lgfortran')
      error('Missing cif2hkl MeX')
    end
  end % ischar
  if isnumeric(this) && length(this) >= 14                   % numeric -> struct
    % [a b c aa bb cc (Z x y z B occ spin charge) .. ], use fixed instrument parameters
    par.cell  = this(1:6); % lattice parameters
    % create the members of the parameter structure
    atom_index=1;
    for j=7:8:length(this)
      vector = this((j+1):min(j+7,length(this)));
      % Z must be an integer value, xyz must be 0<xyz<1
      if this(j) == floor(this(j)) && length(vector) == 8 && all(0<=vector(1:3) && vector(1:3)<=1)
        par.structure.([ atoms(this(j)) num2str(atom_index)]) = vector;
        atom_index = atom_index+1;
      else
        disp([ mfilename ': Unknown trailing numerical arguments. Ignoring' ])
        disp(vector)
      end
    end
    this = par;
  end % is numeric
  if isstruct(this) % struct -> store numeric fields into variable_p, and char into constant_p
    f = fieldnames(this);
    
    for j=1:length(f)
      % check if the name of the field is <atom> optionally followed by a number
      [at,nb] = strtok(f{j}, '0123456789'); % supposed to be an atom, and nb is a 'number' or empty
      %  handle Rietveld parameters
      if any(strcmpi(f{j}, {'Spgr','Spg','Group','SpaceGroup','SubG','SpaceG'}))
        if isnumeric(this.(f{j})), this.(f{j}) = num2str(this.(f{j})); end
        Rietveld.Spgr = this.(f{j});
      elseif any(strncmpi(f{j}, {'struct','atoms'},3))
        Rietveld.structure = this.(f{j});
      elseif any(strncmpi(f{j}, {'CFML_write'},6))
        Rietveld.CFML_write = this.(f{j});
      elseif any(strcmp(at, atoms)) && (isempty(nb) || ~isempty(str2num(nb))) && length(this.(f{j})) >= 3 && length(this.(f{j})) <= 7
        % the name of the field is <atom> optionally followed by a number, and value length is 3-7
        Rietveld.structure.(f{j}) = this.(f{j});
      elseif any(strncmpi(f{j}, {'cell','lattice'},4))
        Rietveld.cell = this.(f{j});
      elseif strncmpi(f{j}, 'title',4)
        Rietveld.title = this.(f{j});
      % handle McStas parameters
      elseif any(strcmpi(f{j},{'ncount','dir','mpi','seed','gravitation','compile','particle','monitors'}))
        mcstas_options.(f{j}) = this.(f{j});
      % handle other parameters (instrument)
      elseif isnumeric(this.(f{j})) && isscalar(this.(f{j}))
        variable_p.(f{j}) = this.(f{j});
      else
        constant_p.(f{j}) = this.(f{j});
      end
    end
  else % end isstruct
    disp([ mfilename ': Invalid input argument of class ' class(this) ':' ])
    disp(this);
  end
end

% check Rietveld structure members =============================================
% Required: cell, Spgr, structure.<atoms>
if ~isstruct(Rietveld)
  error([ mfilename ': No sample parameters defined for Rietveld refinement. See "help rietveld".' ])
end
if ~isfield(Rietveld,'Spgr')
  disp([ mfilename ': Space group not defined. Using default cubic "F m -3 m".' ])
  Rietveld.Spgr = 'F m -3 m'; % default group when not specified
end
if ~isfield(Rietveld, 'cell')
  disp([ mfilename ': lattice constants not defined. Using default a=b=c=2*pi and 90 deg angles.' ])
  Rietveld.cell = [ 2*pi 2*pi 2*pi 90 90 90 ];
end
if ~isfield(Rietveld, 'title')
  Rietveld.title = '';
end
if ~isfield(Rietveld, 'structure')
  error([ mfilename ': No atom positions defined in the cell (Rietveld.structure). See "help rietveld".' ])
end
if ~isstruct(Rietveld.structure)
  error([ mfilename ': The atom positions should be defined as a named structure with [x y z {Biso occ spin charge}] values. See "help rietveld".' ])
end
if ~isfield(Rietveld,'CFML_write')
  Rietveld.CFML_write = 'reflections.laz';
  disp([ mfilename ': Intermediate reflection file not defined. Using default CFML_write=''' Rietveld.CFML_write '''.' ])
end
if isempty(instrument)
  instrument = 'templateDIFF';
end
f = fieldnames(Rietveld.structure);
nb_atoms = length(f);

% build the list of parameters for 'p' and the Guess value =====================
% p(1:6): cell [a b c aa bb cc]
Parameters = {'Rietveld_a','Rietveld_b','Rietveld_c','Rietveld_alpha','Rietveld_beta','Rietveld_gamma'};
Guess      = Rietveld.cell;
Guess(4:6) = mod(Guess(4:6),360);
% p(7:(nb_atoms*7+6)): atoms [x y z {Biso occ spin charge}]
for index=1:length(f)
  this       = Rietveld.structure.(f{index});
  this(1:3)  = min(1,max(0,this(1:3)));
  
  % add optional/missing values {Biso occ spin charge}
  if length(this) < 4, this = [ this 0 ]; end
  if length(this) < 5, this = [ this 1 ]; end
  if length(this) < 6, this = [ this 0 ]; end
  if length(this) < 7, this = [ this 0 ]; end
  this(5) = min(1,max(0,this(5)));
  % setup the Guess for the atom
  Guess      = [ Guess(:) ; this(:) ];
  % add parameter names per atom type
  Parameters{end+1} = [ 'Rietveld_' f{index} '_x' ];
  Parameters{end+1} = [ 'Rietveld_' f{index} '_y' ];
  Parameters{end+1} = [ 'Rietveld_' f{index} '_z' ];
  Parameters{end+1} = [ 'Rietveld_' f{index} '_Biso' ];
  Parameters{end+1} = [ 'Rietveld_' f{index} '_Occ' ];
  Parameters{end+1} = [ 'Rietveld_' f{index} '_Spin' ];
  Parameters{end+1} = [ 'Rietveld_' f{index} '_Charge' ];
end
% add other parameters from 'variable_p'
if isstruct(variable_p)
  f = fieldnames(variable_p);
  for index=1:length(f)
    Parameters{end+1} = f{index};
    Guess = [ Guess variable_p.(f{index}) ];
  end
end
disp([ mfilename ': Assembling Rietveld model ' Rietveld.title ' with ' num2str(length(Guess)) ' parameters, and instrument ' instrument ]);

y.Name       = strtrim([ Rietveld.title ' Rietveld refinement [' mfilename ']' ]);
y.Description= strtrim([ Rietveld.title ' Rietveld refinement using virtual experiment diffractometer with McStas.' ]);
y.Guess      = Guess;
y.Parameters = Parameters;

% assemble the core Expression of the model ====================================
% start to write the CFL file (temporary file name)
Expression = { ...
  'tmp=tempname;', ...
  'fid=fopen(tmp,''w'');', ...
  'fprintf(fid,''! FullProf/CrysFML file format\n'');', ...
  [ 'fprintf(fid,''Title  ' Rietveld.title ' sample\n'');' ], ...
  'fprintf(fid,''!      a           b           c          alpha   beta    gamma\n'');', ...
  'fprintf(fid,''Cell   ''); fprintf(fid,''%f '', p(1:6)); fprintf(fid,''\n'');', ...
  'fprintf(fid,''!     Space Group\n'');', ...
  [ 'fprintf(fid,''Spgr  ' Rietveld.Spgr '\n'');' ], ...
  'fprintf(fid,''!                X        Y       Z     B       occ       Spin  Charge\n'');' ...
  };
% add the Atom list
f = fieldnames(Rietveld.structure);
for index=1:nb_atoms
  i1 = 7+(index-1)*7; % index of <atoms> block in 'p', with 7 values each
  i2 = i1+6;
  lab = upper(strtok(f{index},'0123456789 '));
  Expression{end+1} = [ 'fprintf(fid,''Atom  ' f{index} ' ' lab '''); fprintf(fid,''%g '',p(' ...
    num2str(i1) ':' num2str(i2) ')); fprintf(fid,''\n'');' ];
end
% close the file and request cif2hkl
Expression{end+1} = 'fclose(fid);';

% generate hklF2 (cif2hkl)
Expression{end+1} = [ 'cif2hkl(tmp,''' Rietveld.CFML_write ''');' ];

% delete temporary file
Expression{end+1} = 'delete(tmp);';

% start to build the instrument related Expression
Expression{end+1} = [ 'instrument = ' mat2str(instrument) ';' ];
% additional fit parameters into instrument parameters
if isstruct(variable_p)
  f = fieldnames(variable_p);
  for index=1:length(f)
    i1 = 7+(nb_atoms-1)*7+6+index;  % index of additional variable parameter
    Expression{end+1} = sprintf('instr_pars.%s=p(%i);', f{index}, i1 );
  end
end
% build instrument parameter structure
if isstruct(constant_p)
  f = fieldnames(constant_p);
  for index=1:length(f)
    if isnumeric(constant_p.(f{index})) | ischar(constant_p.(f{index}))
      Expression{end+1} = [ 'instr_pars.' f{index} '=' mat2str(constant_p.(f{index})) ';' ];
    end
  end
end
% build McStas options structure
if isstruct(mcstas_options)
  f = fieldnames(mcstas_options);
  for index=1:length(f)
    if isnumeric(mcstas_options.(f{index})) | ischar(mcstas_options.(f{index}))
      Expression{end+1} = [ 'mcstas_options.' f{index} '=' mat2str(mcstas_options.(f{index})) ';' ];
    end
  end
else
  Expression{end+1} = 'mcstas_options = [];';
end

% evaluate model with mcstas/instrument
Expression{end+1} = [ 'signal = mcstas(instrument, instr_pars, mcstas_options); '];
Expression{end+1} = 'if length(signal) > 1, signal = signal(end); end';

y.Expression = Expression;

% end of core Expression

% create iFunc object
y = iFunc(y);

% check model Dimension. If 0, then start a short McStas simulation to determine the monitor dimensionality
if y.Dimension == 0
  % evaluate 'signal' with Guess parameters
  disp([ mfilename ': Determining the model dimension... (please wait)' ]);
  signal      = feval(y,y.Guess);
  y.Dimension = ndims(signal);
  disp([ mfilename ': Setting dimension of the Rietveld model ' Rietveld.title ' to ' num2str(y.Dimension) ]);
end

