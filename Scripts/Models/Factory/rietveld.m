function y = rietveld(sample, instr, varargin) 
% model=rietveld(sample, instrument, ....) Rietveld refinement of powder/single crystal
%
% This function builds a fit model from:
%   * a sample structure
%   * a McCode instrument model (usually a diffractometer)
% The model can then be used for refinement as a usual fit model:
%   fits(model, data_set, parameters, options, constraints)
%
% MODEL CREATION
%
% Any powder structure can be entered, and there is no limitation on the number
%   of parameters/atoms in the cell. The powder/crystal model takes into account 
%     Biso (thermal motions), charge, occupancy, spin
%   Equally, the instrument model can be of any complexity, including sample
%   environments, other sample phases (incl. amorphous), multi-dimensional PSD 
%   detectors, ... It is possible to specify which monitor file to use from 
%   the McCode instrument. The instrument resolution function is not computed
%   using the legacy Caglioti formalism, but is fully convoluted with the
%   sample component from the McStas instrument description.
%   The default, and recommended, instrument is templateDIFF, but others are possible
%   including TOF-diffractometers, PSD-diffractometers, and even Laue.
%
% Once set, the fit parameters may be constraint with either the usual 'constraints'
%   argument to 'fits', or by setting the model constraints such as:
%     model.Sample_a     = 'fix';     % fix 'a' lattice parameter
%     model.Sample_alpha = [80 110];  % restrict alpha lattice angle in degrees
% 
% Parameters of the model can be entered as any of:
% * a structure with one field 'structure.<atom>' per atom in the cell, named from the atom, e.g.
%     p.structure.Al1=[x y z {Biso occ spin charge}] where {Biso occ spin charge} are optional
%     p.Al1          =[...] can also be used for a compact specification (see exemple below).
%   plus additional fields
%     p.cell=[a b c alpha beta gamma] lattice cell parameters (Angs and deg)
%     p.Spgr='space group' as Number, Hall or Hermman-Mauguin
%
%   In addition, the 'CFML_write' field allows to set the name of the intermediate
%     reflection file written by the CrysFML routines, and usable by the PowderN
%     and Isotropic_Sqw McStas components. The default intermediate file is
%     'reflections.laz'.
%   The optional 'mode' argument to the sample description can be 'p' for powders
%     (default) or 'x' for single crystal.
%
% * a vector array where Z is the atomic number
%     p=[ a b c alpha beta gamma (Z x y z Biso occ spin charge) (...) ... ]
%
% * a CIF/CFL/PCR/shellX file name from which structure is extracted.
%
% * a McCode instrument description file name (.instr extension) to specify
%     the instrument model to use. The default instrument is 'templateDIFF.instr'
%
% * a char with label/value pairs separated by ';' to build a structure.
%
% * the string 'defaults' or 'gui' to use a default model, or pop-up dialogue box.
%
% * All additional arguments are sent to the instrument model.
%
% MODEL EVALUATION
%
% The dimensionality and axes of the model is that of the monitor in the instrument, 
% and its parameters are that from the sample with that of the instrument.
%
% Example: refine a NaCaAlF powder structure with the templateDIFF McStas instrument
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
%    f = rietveld(Sample, 'templateDIFF.instr',' Powder=reflections.laz; lambda="2.36"');
% you may plot the model (using McStas) in the [10 170] angular range:
%    plot(f, f.Guess, linspace(-.2,.2,50),linspace(10,170,300));
% Then 'f' is an iFunc Rietveld model. Then perform the Rietveld refinement:
%    p = fits(measurement, f); % where measurement holds a measured powder diffractogram
%
% See also: mccode, iFunc, iData, iFunc/fits
% (c) E.Farhi, ILL. License: EUPL.

% To initiate the model, one enters parameters as a structure with members.
% During evaluation, parameters are used as a vector matching initiated fields.

y = [];

if nargin == 0
  % No argument -> GUI
  doc(iData,'Models.html#mozTocId344379');  % doc for Rietveld
  NL = sprintf('\n');
  prompt = { [ '{\bf Material structure}' NL ...
    'you should enter a {\color{blue}CIF}, {\color{blue}CFL (FullProf)} or {\color{blue}ShelX} file path, e.g. [ ifitpath ''Data/Na2Ca3Al2F14.cfl'' ].'  ], ...
    [ '{\bf McStas instrument}' NL ...
    'you should enter a {\color{blue}.instr} McStas diffractometer file path, e.g. templateDIFF.instr'  ], ...
    [ '{\bf Additional McStas instrument parameters}' NL ...
    'Any string specifying the {\color{blue}instrument or Mcstas} configuration, such as' NL ...
    '  Powder=reflections.laz; lambda="2.36"; monitor=BananaTheta;' NL ...
    'To use MPI, you may need to force recompilation of the instrument with ' NL ...
    '  ... ; compile=1; mpi=8;' ]
  };
  dlg_title = 'iFit: Model: Rietveld/McStas';
  defAns    = {'[ ifitpath ''Data/Na2Ca3Al2F14.cfl'' ]', 'templateDIFF.instr','Powder=reflections.laz; lambda="2.36"; monitor=BananaTheta; compile=1; mpi=8'};
  num_lines = [ 1 ];
  op.Resize      = 'on';
  op.WindowStyle = 'normal';   
  op.Interpreter = 'tex';
  answer = inputdlg(prompt, dlg_title, num_lines, defAns, op);
  if isempty(answer), 
    return; 
  end
  % extract sample, instrument and more parameters
  try
    sample = eval(answer{1});
  catch
    sample = answer{1};
  end
  try
    instr = eval(answer{2});
  catch
    instr = answer{2};
  end
  try
    options = eval(answer{3});
  catch
    options = answer{3};
  end
  varargin = { options };
  
elseif ischar(sample) && strcmp(sample,'defaults')
  sample = [ ifitpath 'Data/Na2Ca3Al2F14.cfl' ];
  instr  = 'templateDIFF.instr';
  varargin = { 'Powder=reflections.laz; lambda="2.36"; monitor=BananaTheta; mpi=8' };
elseif ischar(sample) && strcmp(sample,'identify')
  y = iFunc;
  y.Name       = [ 'Rietveld refinement [' mfilename ']' ];
  y.Expression = '[]; % dummy code so that it"s not empty p(1)';
  y.Dimension  = -2; % typical but can be something else, e.g. 1-3D
  return
end

if nargin < 2, instr = ''; end

if isempty(instr)
  instr = 'defaults';
  disp([ mfilename ': Instrument description file not defined, using ''' instr '''.' ]);
end


% ==============================================================================
% analyze and check the input structure
CFL       = rietveld_input_atoms(sample, varargin{:});

% ==============================================================================
% build the cif -> laz/lau file model (value is 0)
model_cif = rietveld_cif2hkl(CFL, instr);

% ==============================================================================
% build the McCode model
model_mccode = mccode(instr, varargin{:});

% ==============================================================================
y = model_cif + model_mccode;
disp([ mfilename ': Built model ' y.Name ]);

end % rietveld


% ==============================================================================
function y = rietveld_cif2hkl(CFL, instr)

  f        = fieldnames(CFL.structure);
  nb_atoms = length(f);

  % build the list of parameters for 'p' and the Guess value =====================
  % p(1:6): cell [a b c aa bb cc]
  Parameters = {'Sample_a','Sample_b','Sample_c','Sample_alpha','Sample_beta','Sample_gamma'};
  Guess      = CFL.cell;
  Guess(4:6) = mod(Guess(4:6),360);
  % p(7:(nb_atoms*7+6)): atoms [x y z {Biso occ spin charge}]
  for index=1:length(f)
    this       = CFL.structure.(f{index});
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
    Parameters{end+1} = [ 'Sample_' f{index} '_x' ];
    Parameters{end+1} = [ 'Sample_' f{index} '_y' ];
    Parameters{end+1} = [ 'Sample_' f{index} '_z' ];
    Parameters{end+1} = [ 'Sample_' f{index} '_Biso' ];
    Parameters{end+1} = [ 'Sample_' f{index} '_Occ' ];
    Parameters{end+1} = [ 'Sample_' f{index} '_Spin' ];
    Parameters{end+1} = [ 'Sample_' f{index} '_Charge' ];
  end
  
  disp([ mfilename ': Assembling Rietveld model ' CFL.title ' with ' num2str(length(Guess)) ' parameters, and instrument ' instr ]);

  % assemble the core Expression of the model ====================================
  % start to write the CFL file (temporary file name)
  Expression = { ...
    'tmp = [ tempname ''.cfl'' ];', ...
    'fid=fopen(tmp,''w'');', ...
    'fprintf(fid,''! FullProf/CrysFML file format\n'');', ...
    [ 'fprintf(fid,''Title  ' CFL.title ' sample\n'');' ], ...
    'fprintf(fid,''!      a           b           c          alpha   beta    gamma\n'');', ...
    'fprintf(fid,''Cell   ''); fprintf(fid,''%f '', p(1:6)); fprintf(fid,''\n'');', ...
    'fprintf(fid,''!     Space Group\n'');', ...
    [ 'fprintf(fid,''Spgr  ' CFL.Spgr '\n'');' ], ...
    'fprintf(fid,''!                X        Y       Z     B       Occ       Spin  Charge\n'');' ...
    };
  % add the Atom list
  f = fieldnames(CFL.structure);
  for index=1:nb_atoms
    i1 = 7+(index-1)*7; % index of <atoms> block in 'p', with 7 values each
    i2 = i1+6;
    lab = upper(strtok(f{index},'0123456789 '));
    Expression{end+1} = [ 'fprintf(fid,''Atom  ' f{index} ' ' lab ' ''); fprintf(fid,''%g '',p(' ...
      num2str(i1) ':' num2str(i2) ')); fprintf(fid,''\n'');' ];
  end
  % close the file and request cif2hkl
  Expression{end+1} = 'fclose(fid);';

  % generate hklF2 (cif2hkl)
  if isfield(CFL, 'mode')
    Expression{end+1} = [ 'cif2hkl(tmp,''' CFL.CFML_write ''',[],''' CFL.mode ''');' ];
  else
    Expression{end+1} = [ 'cif2hkl(tmp,''' CFL.CFML_write ''');' ];
  end

  % delete temporary file
  Expression{end+1} = 'delete(tmp);';
  Expression{end+1} = 'signal = 0;';  % to be able to add with the instrument

  y.Name       = strtrim([ CFL.title ' Rietveld refinement [' mfilename ']' ]);
  y.Description= strtrim([ 'Rietveld refinement using McCode virtual experiment ' instr ]);
  y.Guess      = Guess;
  y.Parameters = Parameters;
  y.ParameterValues = Guess;
  y.Expression = Expression;

  % create iFunc object
  y = iFunc(y);
  
  % store cell and reciprocal basis 'B'
  y.UserData.cell = Guess;
  alpha=y.UserData.cell(4);
  beta =y.UserData.cell(5);
  gamma=y.UserData.cell(6);
  a_vec=y.UserData.cell(1)*[1; 0; 0];
  b_vec=y.UserData.cell(2)*[cosd(gamma); sind(gamma); 0];
  c1=cosd(beta);
  c2=(cosd(alpha)-cosd(gamma)*cosd(beta))/sind(gamma);
  c3=sqrt(1-c1^2-c2^2);
  c_vec=y.UserData.cell(3)*[c1; c2; c3;];
  V=dot(a_vec,cross(b_vec,c_vec));
  % reciprocal basis, as columns
  y.UserData.reciprocal_cell=2*pi*[cross(b_vec,c_vec) cross(c_vec,a_vec) cross(a_vec,b_vec)]/V; 
  y.UserData.volume = V;

  % restraint structure parameters: xyz in [0 1], angles in [0 180]
  %                                 Biso in [0 10], Occ in [0 1]
  %                                 Spin Charge in [-10 10]
  xlim(y, regexp(y.Parameters, '_x\>|_y\>|_z\>|_Occ\>'), [0 1]); % \> = end of Parameter names
  xlim(y, {'Sample_alpha','Sample_beta','Sample_gamma'}, [0 180]);
  xlim(y, regexp(y.Parameters, '_Biso\>'), [0 10]);
  xlim(y, regexp(y.Parameters, '_Spin\>|_Charge\>'), [-10 10]);

end % rietveld_cif2hkl


function CFL = rietveld_input_atoms(varargin)

  CFL            = [];

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
    if isempty(this), continue; end
    if ischar(this) && ~isempty(dir(this))  % char -> file -> struct OR instrument
      % call cif2hkl
      this = read_cif(this);
      if ~isempty(this), CFL = this; end
    end % ischar
    if ischar(this)   % char -> struct
      par = str2struct(this); 
      if isstruct(par), this = par; end
    end
    % now 'this' should not be a char anymore (has been transformed to CFL, numeric or struct)
    if ischar(this)
      disp([ mfilename ': Unknown char argument ' this '. Ignoring' ])
      continue
    end
    
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
    if isstruct(this) % struct -> may be a CFL given as a structure ?
      [CFLnew, this] = read_cif(this);
      if isempty(CFL) && ~isempty(CFLnew), CFL = CFLnew; end
    else % end isstruct
      disp([ mfilename ': Invalid input argument of class ' class(this) ':' ])
      disp(this);
    end
  end % for
  
  % check CFL structure members =============================================
  % Required: cell, Spgr, structure.<atoms>
  if ~isstruct(CFL)
    error([ mfilename ': No sample parameters defined for Rietveld refinement. See "help rietveld".' ])
  end
  if ~isfield(CFL, 'cell')
    disp([ mfilename ': lattice constants not defined. Using default a=b=c=2*pi and 90 deg angles.' ])
    CFL.cell = [ 2*pi 2*pi 2*pi 90 90 90 ];
  end
  if ~isfield(CFL, 'title')
    CFL.title = '';
  end
  if ~isfield(CFL, 'structure')
    error([ mfilename ': No atom positions defined in the cell (CFL.structure). See "help rietveld".' ])
  end
  if ~isstruct(CFL.structure)
    error([ mfilename ': The atom positions should be defined as a named structure with [x y z {Biso occ spin charge}] values. See "help rietveld".' ])
  end
  if ~isfield(CFL,'CFML_write')
    CFL.CFML_write = 'reflections.laz';
    disp([ mfilename ': Intermediate reflection file not defined. Using default CFML_write=''' CFL.CFML_write '''.' ])
  end
  if ~isfield(CFL,'Spgr')
    disp([ mfilename ': Space group not defined. Using default cubic "F m -3 m".' ])
    CFL.Spgr = 'F m -3 m'; % default group when not specified
  end
end % rietveld_input_atoms

