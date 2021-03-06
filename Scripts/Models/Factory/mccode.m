function y = mccode(instr, options, parameters)
% y = mccode(instr, options, parameters) : McCode (McStas/McXtrace) instrument
%
%   iFunc/mccode a McCode instrument
%     y=model instrument
%
% MODEL CREATION:
% ------------------------------------------------------------------------------
% mccode(description)
%       creates a model with specified McCode instrument
%       The instrument may be given as an '.instr' McCode description, or directly
%       as an executable. Distant URL (ftp, http, https) can also be used.
% mccode('')
%       requests a McCode file (*.instr,*.out) with a file selector.
% mccode('gui')   and 'mccode'  alone.
%       list all available instruments in a list for a selection.
% mccode('defaults')
%       uses templateDIFF.instr neutron powder diffractometer as example.
% mccode('check')
% mccode('check', options)
%       check for the McStas/McXtrace installation. You can indicate explicitly
%       the location of the 'mcrun/mxrun'  command with e.g.
%       options.mccode='/usr/local/bin/mcrun'. The 'mcstas' and 'mcrun/mxrun'
%       commands should be available.
% mccode(description, options) also specifies additional McCode options, e.g.
%   options.dir:         directory where to store results, or set automatically (string)
%                          the last simulation files are stored therein 'sim'.
%                        dir can also be 'pwd' for current, or 'tmp' for temporary location.
%   options.ncount:      number of neutron events per iteration, e.g. 1e6 (double)
%   options.mpi:         number of processors/cores to use with MPI on localhost (integer)
%                          when MPI is available, and mpi options is not given,
%                          all cores are then used.
%   options.machines:    filename containing the list of machines/nodes to use (string)
%   options.seed:        random number seed to use for each iteration (double)
%   options.gravitation: 0 or 1 to set gravitation handling in neutron propagation (boolean)
%   options.monitor:     a single monitor name to read, or left empty for the last (string).
%                        this can be a wildcard expression.
%   options.mccode:      set the executable path to 'mcrun' (default, neutrons) or 'mxrun' (xrays)
%   options.mpirun:      set the executable path to 'mpirun'. You can set
%                        it to 'none' to not use MPI.
%   options.compile:     0 or 1 to compilation the executable. The default is to compile.
%
%   All options are stored and assignable in model.UserData.options.
%   options can also be given as a string, e.g. 'ncount=1e6; monitor=*Theta*; compile=1'
%   the 'monitor' option can also include further expressions, such as:
%     options.monitor='*Theta*; signal=max(signal)/std(signal)^2;'
%
%   When the instrument file name contains a path specification, it is used for
%   compilation and execution (recommended). If it is used as a single file name,
%   or the target directory as no read, write or execute permissions,
%   a temporary directory is created and the instrument is copied there, which may
%   cause issues if local components are to be used.
%
% mccode(description, options, parameters)
%   Specifies the instrument parameters values to use as default. These values can
%   be given as a string e.g. 'QM=1; lambda=2.36' or a structure.
%
% The instrument parameters of type 'double' are used as model parameters. Other
% parameters (e.g. of type string and int) are stored in UserData.Parameters_Constant
% You may still specify their value when using e.g. syntax with a string or structure:
%   feval(model, 'lambda=2.36; Powder="Al.laz"')
%   iData(model, 'lambda=2.36; Powder="Al.laz"')
%   iData(model, struct('lambda',2.36,'Powder','Al.laz'))
%
% The options ncount, seed, gravitation, monitor can be changed for the model
% evaluation, with e.g.:
%   model.UserData.options.ncount     =1e5;
%   model.UserData.options.gravitation=1;
%   model.UserData.options.monitor    ='*Theta*';
%
% Additional information is stored in the model.UserData, such as the instrument
% source, which you may view with:
%   TextEdit(model.UserData.instrument_source)
%
% Example: model = mccode('templateDIFF');
%
% MODEL EVALUATION:
% ------------------------------------------------------------------------------
% model(p)
%   evaluates the model with given parameters (vector, cell, structure). Only
%   scalar/double parameters of the instrument can be varied. Other parameters are kept fixed.
% model(p, nan)
%   evaluates the model and return the raw McCode data set (monitor).
% model(p, x,y,...)
%   evaluates the model and interpolates the McCode data set onto given axes.
%
% The model 'value' is the last monitor, or set from
%   model.UserData.options.monitor
% It can be converted to an iData with iData(model, ...)
%
% The raw monitors of the last simulation are stored as iData objects in
%   model.UserData.monitors
% You can plot them all with e.g.:
%   subplot(model)
%
% input:  p: variable instrument parameters (double, struct, char)
%            p = [ double_type_instrument_parameters ]
%         x,y,...: axes (double)
%
% output: y: monitor value
% ex:     model =mccode('templateDIFF');
%         signal=iData(model, [], linspace(-10,100,100));
%         signal=iData(model, [], nan); % to get the raw monitor
%         signal=iData(model, 'lambda=2.36; Powder=Al.laz', nan); % to get the raw monitor
%
% MODEL GEOMETRY
% ------------------------------------------------------------------------------
% To view the model geometry, use the plot routine:
%   model = mccode('templateDIFF');
%   plot(model)
%   plot(model, parameters)
%
% MODEL PARAMETER SCAN
% ------------------------------------------------------------------------------
% It is possible to scan model parameters when using vectors as value for the
% parameters. To achieve that, the parameter values must be given as a named
% structure.
% For instance, to scan the RV parameter in the templateDIF instrument, use:
%   model = mccode('templateDIFF');
%   p.RV= [ 0:.25:2 ];      % from 0 to 2 by steps of .25
%   v=iData(model, p, nan); % we want the raw monitors as iData sets.
%   subplot(v);
%
% MODEL OPTIMISATION
% ------------------------------------------------------------------------------
% To optimise instrument parameters, you should first fix the non-varying
% parameters, and possibly bound the others. Then the optimiser is launched
% with e.g. 'fmax' or 'fmin':
%   model = mccode('templateDIFF');  % maximize
%   fix(model, 'all'); model.RV='free';
%   model.RV=[0 1 2];        % bounds and starting value
%   [p,c,msg,output] = fmax( model, [], 'OutputFcn=fminplot', nan)  % return the optimal parameters using the raw monitors
%
% Version: $Date$ $Version$ $Author$
% See also iFunc, iFunc/fits, iFunc/plot, iFunc/feval, mcstas
%          <a href="http://www.mcstas.org">McStas</a>, <a href="http://www.mccode.org">McCode</a>
%

% check for McCode executable
persistent mccode_present

y = [];
% get options
if nargin == 0
  instr = '';
end

if nargin > 1
  if ischar(options)
    options = str2struct(options);
  elseif ~isstruct(options)
    options=[];
  end
else
  options=[];
end
if isempty(options)
  options = struct();
end
options = instrument_parse_options(options);

if nargin < 3,
  parameters = [];
end
if ischar(parameters)
  parameters = str2struct(parameters);
end

% check for McCode, only the fisrt time
if isempty(mccode_present) || (nargin>0 && strcmp(instr, 'check'))
  mccode_present = mccode_check(options);
end

if strcmp(instr, 'check') || strcmp(instr, 'config')
  y = mccode_present;
  return
end
% get the instrument. If not available, search for one in a McCode installation
% stop if not found

% when nothing given, ask user. a list selector with all found instruments.
if ischar(instr) && strcmp(instr, 'gui')
  instr= mccode_search_instrument('.instr', options.dir);
  [selection] = listdlg('PromptString', ...
    {'Here are all found instruments on your system.'; ...
     'Select a McStas/McXtrace instrument to load'},...
    'Name','McCode: Select a McStas/McXtrace instrument to load', ...
    'SelectionMode','single',...
    'ListSize', [500 300], ...
    'ListString',instr);
  if isempty(selection),    return; end
  options.instrument = instr{selection};
elseif ischar(instr) && strcmp(instr,'identify')
  y = iFunc;
  y.Name       = [ 'McCode Monte-Carlo neutron/X-ray instrument [' mfilename ']' ];
  y.Expression = '[]; % dummy code so that it"s not empty p(1)';
  y.Dimension  = -2; % typical but can be something else, e.g. 1-3D
  return
else
  % empty choice: pop-up a file selector
  if isempty(instr),
    filterspec = {'*.*', 'All files (*.*)' ; ...
                  '*.instr', 'McCode instrument (*.instr)' ; ...
                  '*.out;*.exe','McCode compiled instrument executable (*.out, *.exe)' };
    [filename, pathname] = uigetfile(filterspec,'Select a McStas/McXtrace instrument to load');
    if isempty(filename),    return; end
    if isequal(filename, 0), return; end
    options.instrument = fullfile(pathname, filename);
  else
    if strcmp(instr, 'defaults')
      instr='templateDIFF.instr';
      if ~isfield(options,'dir') options.dir='pwd'; end
    end
%    elseif strcmp(instr, 'identify')
%      y = mccode('defaults');
%      y.Name = [ 'McCode instrument [' mfilename ']' ];
%      return;
%    end
    options.instrument = mccode_search_instrument(instr, options.dir);
  end
end

if isempty(options.instrument), return; end
if iscell(options.instrument), options.instrument = options.instrument{1}; end
disp([ mfilename ': Using instrument: ' options.instrument ] );

% target directory aliases
if strcmp(options.dir, 'pwd')
  options.dir = pwd;
elseif any(strcmp(options.dir, {'tmp','temp','tempname'}))
  options.dir = tempname;
end

% determine if we use a fully qualified path in instrument, then sets it.
[p,f,e] = fileparts( options.instrument );
if ~isempty(p) && isempty(options.dir)
  options.dir = p;
end

% test for directory permissions
[~,attr] = fileattrib(options.dir);
if isstruct(attr) && ~isempty(options.dir) && (~attr.UserRead || ~attr.UserWrite || ~attr.UserExecute)
  disp([ mfilename ': target directory ' options.dir 'is not accessible. Using temporary.' ]);
  options.dir = ''; % will use temporary dir
end
% attempt to create a file there
if ~isempty(options.dir)
  testDir = ['deleteMe_', num2str(floor(rand*1e12))];
  [isWritable,message,messageid] = mkdir(options.dir, testDir);
  if ~isWritable, options.dir = ''; % will use TMPDIR
  else rmdir(fullfile(options.dir, testDir),'s');
  end
end

% use temporary directory to build/assemble parts.
if isempty(options.dir),
  options.dir = tempname;
end
if ~isdir(options.dir)
  mkdir(options.dir);
end

% copy file locally
try
  copyfile(options.instrument, options.dir);
end
options.instrument = fullfile(options.dir,[f e]);

if isfield(options, 'mpi') && options.mpi <= 1
  options.mpirun = ''; % no MPI needed
end
if isempty(options.mpirun)
  options.mpirun = mccode_present.mpirun; % get MPIRUN from system
else
  if strcmp(options.mpirun, 'none')
      options.mpirun=''; % no MPI use
  else
    mccode_present.mpirun = options.mpirun;
  end
end
if isempty(options.mccode)
  options.mccode = mccode_present.mccode;
else
  mccode_present.mccode = options.mccode;
end
% MPI: default to all available core (read from Java)
if ~isempty(options.mpirun) && ~isfield(options,'mpi')
  if usejava('jvm')
    r=java.lang.Runtime.getRuntime;
    % mem_avail   = r.freeMemory;
    options.mpi = r.availableProcessors;
  else try; options.mpi = feature('NumCores'); end; end
end

% test if the given/found instrument is executable. Return empty if source code.
[info,exe] = instrument_get_info(options.instrument); % try --info and --help

% compile the instrument (when McCode present), possibly with MPI from options
% if no McCode executable, the user must provide an executable.
if (isempty(info) && ~isempty(options.mccode)) || options.compile
  % compile. Stop when fails.
  result = instrument_compile(options);
  % identify and test the executable
  [info, exe] = instrument_get_info(options.instrument);
  if isempty(info)
    disp([ mfilename ': ERROR: Compiled instrument ' exe ' from ' options.instrument ' is not executable.'])
    return
  end
  % store the source code, and binary (executable) in UserData
  UserData.instrument_source     = fileread(options.instrument);
elseif ~isempty(info)
  % store the executable, no source code given
  UserData.instrument_source     = '';
else
  disp([ mfilename ': WARNING: No information could be retrieved from ' options.instrument ]);
end
[fid, message] = fileattrib(exe);
if     ~fid || ~isstruct(message), fid=-1;
elseif isfield(message,'UserExecute') && ~message.UserExecute, fid=-1; end
if fid == 1
  [fid, message] = fopen(exe, 'r');
end
if fid == -1
  disp(message);
  disp(info);
  error([ mfilename ': the executable ' exe ' is not accessbile. Perhaps you have a system command with the same name ? Change instrument name ?' ]);
end
UserData.instrument_executable = uint8(fread(fid, Inf)); fclose(fid);
UserData.instrument_exe = exe;

[p,f,e] = fileparts(options.instrument); options.instrument = [ f e ];
[p,f,e] = fileparts(UserData.instrument_exe);        UserData.instrument_exe = [ f e ];

% get default parameters (and test executable). Store them in UserData
UserData.instrument_info    = info;
UserData.options            = options;

% identify fixed parameters (non scalar). Store them in UserData
% we search the 'Param_*' fields in info.
[UserData.Parameters, UserData.Parameters_Constant] = instrument_get_parameters(info);

y.Parameters = fieldnames(UserData.Parameters);
% set starting configuration (Guess)
c = struct2cell(UserData.Parameters);
f = fieldnames(UserData.Parameters);

y.Guess = [];
for index=1:numel(c)
  this = c{index};
  % override defined value in instrument with that given in 'parameters' input arg.
  if isfield(parameters, f{index})
    this = parameters.(f{index});
  end
  if ~isempty(this) && isfinite(this),
       y.Guess(end+1) = this;
  else y.Guess(end+1) = nan; end
end

disp([ mfilename ': Assembling McCode model ' options.instrument ' with ' num2str(length(y.Guess)) ' parameters.'  ]);
if ~isempty(UserData.Parameters_Constant) && ~isempty(fieldnames(UserData.Parameters_Constant))
  disp('  Additional parameters are stored in UserData.Parameters_Constant:')
  disp(UserData.Parameters_Constant)
end

y.Name       = strtrim([ options.instrument ' McCode [' mfilename ']' ]);
NL = sprintf('\n');
y.Description= [ 'McCode virtual experiment ' options.instrument ...
    NL '  Set UserData.options.monitor to specify a given monitor file pattern, ' ...
    NL '  or [] to get the last.' ...
    NL '  Monitors are stored in UserData.monitors' ];
y.UserData   = UserData;

% assemble the Expression:
%     build target directory (put back executable if missing)
%     assemble all parameters and options (command line). Handle MPI and machine list
%     execute
%     get back results (filter monitor)
%     set signal (interp to Axes)
ax = ',x,y,z,t';
y.Expression = { ...
'  UD = this.UserData; options=UD.options;', ...
'  % check directory for execution', ...
'  if ~isempty(options.dir) && isempty(dir(options.dir))', ...
'    try; mkdir(options.dir); end', ...
'  end', ...
'  if isempty(options.dir) || isempty(dir(options.dir))', ...
'    options.dir=tempname;', ...
'    mkdir(options.dir)', ...
'    options.use_tmpdir = true;', ...
'  else options.use_tmpdir = false;', ...
'  end', ...
'  % check if executable is present. Write it if missing.', ...
'  if isempty(dir(fullfile(options.dir, UD.instrument_exe)))', ...
'    fid=fopen(fullfile(options.dir, UD.instrument_exe), ''w'');', ...
'    if fid==-1, error([ ''model '' this.Name '' '' this.Tag '' could not write '' fullfile(options.dir, UD.instrument_exe) ]); end', ...
'    fwrite(fid, UD.instrument_executable);', ...
'    fclose(fid); fileattrib(fullfile(options.dir, UD.instrument_exe),''+x'');', ...
'  end', ...
'  % check if source is present. Write it if missing.', ...
'  if ~isempty(UD.instrument_source) && isempty(dir(fullfile(options.dir, options.instrument)))', ...
'    fid=fopen(fullfile(options.dir, options.instrument), ''w'');', ...
'    if fid==-1, disp([ ''WARNING: model '' this.Name '' '' this.Tag '' could not write '' fullfile(options.dir, options.instrument) ]); end', ...
'    fprintf(fid, ''%s\n'', UD.instrument_source);', ...
'    fclose(fid);', ...
'  end', ...
'  % handle mpi', ...
'  if isfield(options,''mpi'') && ~isempty(options.mpi) && ...', ...
'    ~isempty(options.mpirun) && isnumeric(options.mpi) && options.mpi > 1 && ...', ...
'    (isempty(options.trace) || ~isscalar(options.trace) || ~options.trace)', ...
'    cmd = [ options.mpirun '' -n '' num2str(options.mpi) '' '' ];', ...
'    if isfield(options, ''machines'') && ~isempty(options.machines) && ischar(options.machines)', ...
'      cmd = [ cmd '' -machinefile '' options.machines '' '' ];', ...
'    end', ...
'  else cmd = ''''; options.mpi=1; end', ...
'  % assemble command line', ...
'  if ~isempty(dir(fullfile(options.dir,''sim'')))', ...
'    try; rmdir(fullfile(options.dir,''sim''),''s''); catch ME; disp(getReport(ME)); end; end', ...
'  cmd = [ cmd fullfile(options.dir, UD.instrument_exe) '' --ncount='' num2str(options.ncount) ...', ...
'    '' --dir='' fullfile(options.dir,''sim'') '' '' options.raw ];', ...
'  % handle seed gravitation', ...
'  if isfield(options,''seed'') && ~isempty(options.seed) && isscalar(options.seed)', ...
'    cmd = [ cmd '' --seed '' options.seed ];', ...
'  end', ...
'  if isfield(options,''gravitation'') && ~isempty(options.gravitation) && isscalar(options.gravitation) && options.gravitation', ...
'    cmd = [ cmd '' --gravitation '' ];', ...
'  end', ...
'  if ~isempty(options.trace) && isscalar(options.trace) && options.trace', ...
'    cmd = [ cmd '' --trace '' ];', ...
'  end', ...
'  % add parameters', ...
'  f = this.Parameters;', ...
'  for index=1:numel(f) % variable instrument parameters', ...
'    cmd = [ cmd '' '' f{index} ''='' num2str(p(index)) ];', ...
'  end', ...
'  f = fieldnames(UD.Parameters_Constant);', ...
'  for index=1:numel(f) % constant (other) instrument parameters', ...
'    p0 = UD.Parameters_Constant.(f{index});', ...
'    if isempty(p), continue; end', ...
'    if isnumeric(p0)', ...
'      cmd = [ cmd '' '' f{index} ''='' num2str(p0) ];', ...
'    elseif ischar(p0)', ...
'      cmd = [ cmd '' '' f{index} ''="'' p0 ''"'' ];', ...
'    end', ...
'  end', ...
'  % execute', ...
'  disp(cmd);', ...
'  [status, result] = system(cmd);', ...
'  disp(result)', ...
'  if status', ...
'    error([ ''Model '' this.Name '' '' this.Tag '' failed to execute '' cmd ]);', ...
'  end', ...
'  % read results', ...
'  signal = [];', ...
'  if ~isempty(options.monitor) && ischar(options.monitor)', ...
'    [pattern,expr] = strtok(options.monitor, '' ;:='');', ...
'    signal = iData(fullfile(options.dir,''sim'',pattern));', ...
'    if ~isempty(expr) && any(expr(1)=='' :;='') expr=expr(2:end); end', ...
'    try; eval([ expr '';'' ]); catch; disp([ ''ERROR: '' this.Tag '': mccode: Ignoring faulty monitor expression '' expr]); end', ...
'  end', ...
'  if all(isempty(signal)) || isempty(options.monitor)', ...
'    if ~isempty(dir(fullfile(options.dir,''sim'',''mccode.sim'')))' ...
'      signal = iData(fullfile(options.dir,''sim'',''mccode.sim''));', ...
'    elseif ~isempty(dir(fullfile(options.dir,''sim'',''mcstas.sim'')))' ...
'      signal = iData(fullfile(options.dir,''sim'',''mcstas.sim''));', ...
'    end' ...
'  end', ...
'  if options.use_tmpdir, rmdir(options.dir,''s''); end', ...
'  this.UserData.monitors = signal;', ...
'  if ~isempty(options.monitor) && isnumeric(options.monitor) && numel(signal) > 1', ...
'    signal=signal(options.monitor);', ...
'  end', ...
'  if isa(signal, ''iData'')', ...
'    if numel(signal) > 1, signal = signal(end); end', ...
'    this.Dimension = ndims(signal);', ...
['    ax=''' ax(2:end) ''';' ], ...
'    nd=exist(''t'')+exist(''z'')+exist(''y'')+exist(''x'');' ...
'    if min(nd,this.Dimension)>0, ax=eval([ ''{'' ax(1:(2*min(nd,this.Dimension))) ''}'']); else ax={}; end', ...
'    if ~isempty(ax) && exist(''x'') && ~isempty(x) && ~all(isnan(x(:))), signal = interp(signal, ax{:});', ...
'    else x=getaxis(signal,1); y=getaxis(signal,2); z=getaxis(signal,3); t=getaxis(signal,4); end;', ...
'  end', ...
'' ...
};

y.ParameterValues = y.Guess;

% create the iFunc model
y = iFunc(y);
% make it an iFunc sub-class
y = iFunc_McCode(y);

% evaluate the model to get its monitors, and derive the dimensionality
% all parameters must be defined to get execute with default values
%
% at this stage, 'signal' contains all monitors
disp([ mfilename ': Determining the model dimension...' ]);
ncount = options.ncount;
y.UserData.options.ncount = 1e2;  % fast as we only need the size and position
y.UserData.monitors    = [];

% evaluate the model. This may fail if some ERROR is raised by the
% instrument itself (e.g. missing input parameter, wrong configuration,
% ...)
try
    [signal,y] = feval(y,y.Guess, nan);
    signal = y.UserData.monitors;     % get all monitors
catch
    signal = [];
end
y.UserData.options.ncount = ncount; % restore initial ncount
if isa(signal,'iData')
  if numel(signal) > 1  % get the last one
    signal=signal(end);
  end
  y.Dimension = ndims(signal);

  if numel(y.UserData.monitors) > 1
    y.Description = [ y.Description sprintf('\n  Available monitors:') ];
    for index=1:numel(y.UserData.monitors)
      y.Description = [ y.Description ...
        sprintf('\n    * %s', get(y.UserData.monitors(index),'Component')) ];
    end
  end
  for index=1:y.Dimension
    x = getaxis(signal, index); x=x(:);
    y.Description = [ y.Description ...
      sprintf('\n') char(signal) ...
      sprintf('\n  Axis %i "%s" label is "%s", range [%g:%g]', index, ax(2*index), label(signal, index), min(x), max(x)) ];
  end
elseif isnumeric(signal)
  if isscalar(signal),     y.Dimension = 0;
  elseif isvector(signal), y.Dimension = 1;
  else y.Dimension = ndims(signal); end
end

disp(y.Description);
if ~isdeployed
  disp([ mfilename ': built model ' y.Name ' in <a href="' y.UserData.options.dir '">' y.UserData.options.dir '</a>' ])
else
  disp([ mfilename ': built model ' y.Name ' in ' y.UserData.options.dir ])
end


% ==============================================================================
function options = instrument_parse_options(options)
  if ~isfield(options,'instrument'), options.instrument = ''; end
  if ~isfield(options,'ncount'),     options.ncount     = 1e6; end
  if ~isfield(options,'dir'),        options.dir        = ''; end
  if ~isfield(options,'monitor'),    options.monitor    = ''; end
  if ~isfield(options,'gravitation'),options.gravitation= 0; end
  if ~isfield(options,'seed'),       options.seed       = []; end
  if ~isfield(options,'mccode'),     options.mccode     = []; end
  if ~isfield(options,'mpirun'),     options.mpirun     = 'mpirun'; end
  if ~isfield(options,'trace'),      options.trace      = ''; end
  if ~isfield(options,'compile'),    options.compile    = 1; end
  if ~isfield(options,'raw'),        options.raw        = ''; end

% ------------------------------------------------------------------------------
function present = mccode_check(options)
% check if McCode (mcstas or mcxtrace) is present

  % required to avoid Matlab to use its own libraries
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ; DISPLAY= ; ';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ;  DISPLAY= ; ';
  else           precmd = ''; end

  present.mccode = '';
  for totest = { options.mccode, 'mcrun-pl','mxrun-pl','mcrun','mxrun' }
    if ~isempty(present.mccode), break; end
    for ext={'','.pl','.py','.bat','.exe','.out'}
      % look for executable and test with various extensions
      [status, result] = system([ precmd totest{1} ext{1} ]); % usually 127 indicates 'command not found'
      if any(status == [ -1 1 255]) && (~ispc || isempty(strfind(result, [ '''' totest{1} ext{1} '''' ])))
        present.mccode = [ totest{1} ext{1} ];
        break
      end
    end
  end

  if isempty(present.mccode)
    disp([ mfilename ': WARNING: ' options.mccode ' McStas/McXtrace executable is not installed. Get it at www.mccode.org' ]);
    disp('  The model can still be created if the instrument is given as an executable.')
    disp('You may try "mccode check" after extending the PATH to the MCSTAS executable location with e.g.')
    disp('  setenv(''PATH'', [getenv(''PATH'') '':/usr/local/bin'' '':/usr/bin'' '':/usr/share/bin'' ]);');
    disp('You may as well indicate the McCode executable with e.g. options.mccode=''/usr/local/bin/mccrun''.' );
    disp('It may also be required so set the MCSTAS_CC environment variable to an active C compiler.')
  else
    disp([ '  McCode          (http://www.mccode.org) as "' present.mccode '"' ]);
  end

  % test for mpirun
  present.mpirun = '';
  for calc={options.mpirun, 'mpirun', 'mpiexec'}
    % now test executable
    [st,result]=system([ precmd 'echo "0" | ' calc{1} ]);
    if any(st == [ -1 0:2 ]) && (~ispc || isempty(strfind(result, [ '''' calc{1} '''' ])))
        present.mpirun=calc{1};
        st = 0;
        disp([ '  MPI             (http://www.openmpi.org) as "' present.mpirun '"' ]);
        break;
    end
  end

  % test for C compiler
  present.cc = '';
  for try_cc={getenv('CC'),'cc','gcc','ifc','pgcc','clang','tcc'}
    if ~isempty(try_cc{1})
      [status, result] = system([ precmd try_cc{1} ]);
      if status == 4 || ~isempty(strfind(result,'no input file'))
        present.cc = try_cc{1};
        disp([ '  CC              C compiler as "' present.cc '"' ]);
        break;
      end
    end
  end
  if isempty(present.cc)
    disp([ mfilename ': ERROR: C compiler is not available from PATH:' ])
    disp(getenv('PATH'))
    disp([ mfilename ': You may have to extend the PATH with e.g.' ])
    disp('setenv(''PATH'', [getenv(''PATH'') '':/usr/local/bin'' '':/usr/bin'' '':/usr/share/bin'' ]);');

    warning('%s: Can''t find a valid C compiler. Install any of: gcc, ifc, pgcc, clang, tcc\n', ...
    mfilename);
  end

% ------------------------------------------------------------------------------
function instr = mccode_search_instrument(instr, d)

  % get/search instrument
  % check if the instrument exists, else attempt to find it

  if strncmp(instr,'http://',7) || strncmp(instr,'https://',8) || strncmp(instr,'ftp://',6)
    tmpfile = tempname;
    % Keep file extension, may be useful for iData load
    [filepath,name,ext] = fileparts(instr);
    tmpfile = [tmpfile ext];
    use_wget = false;
    if ~usejava('jvm')
      use_wget = true;
    else
      % access the net. Proxy settings must be set (if any).
      try
        % write to temporary file
        tmpfile = urlwrite(instr, tmpfile);
      catch ME
        use_wget = true;
      end
    end
    if use_wget
      % Fall back to using wget
      cmd = ['wget ' instr ' -O ' tmpfile]; disp(cmd)
      [status, result] = system(cmd);
      if status
        disp(result);
        error([ mfilename ': Can not get URL ' instr ]);
      end
    end
    instr = tmpfile;
  end


  if ~isempty(instr)
    index = dir(instr);
  else return;
  end

  if ~isempty(index), return; end % given file is fully qualified

  for ext={'.instr','.out','.exe',''}
    out = [ instr ext{1} ];
    % check for instrument in McStas/McXtrace libraries
    search_dir = { '.', d, getenv('MCSTAS'), getenv('MCXTRACE'), ...
      '/usr/local/lib/mc*', 'C:\mc*', '/usr/share/mcstas/'};
    if isempty(index)
      % search the instrument recursively in all existing directories in this list
      index = getAllFiles(search_dir, out);
      % check if we have more than one match
      if ~isempty(index)
        if numel(index) > 1
          % filter with extension
          found = {};
          for i=1:numel(index)
            [p,f,e] = fileparts(index{i});
            if strcmp(e, '.instr'), found{end+1} = index{i}; end
          end
          if ~isempty(found), index=found; end
          disp([ mfilename ': Found instruments:' ] );
          fprintf(1,'  %s\n', index{:});
        end
        instr = index;
        return
      end
    end
  end

  if isempty(index)
    disp([ mfilename ': ERROR: Can not find instrument ' instr ]);
    instr = '';
  end

% ------------------------------------------------------------------------------
% function to search for a file recursively
function fileList = getAllFiles(dirName, File, recursive)

  if nargin < 3, recursive = true; end

  % allow search in many directories
  if iscell(dirName)
    fileList = {};
    for d=1:length(dirName)
      if strcmp(dirName{d}, '.')
        fileList1=getAllFiles(dirName{d}, File, false); % not recursive for '.' as it may be long
      else
        fileList1=getAllFiles(dirName{d}, File);
      end
      if ~isempty(fileList1)
        if ischar(fileList1)
          fileList{end+1} = fileList1;
        else
          fileList = { fileList{:} fileList1{:} };
        end
      end
    end
    return
  end

  dirData = dir(dirName);                 % Get the data for the current directory
  fileList= {};
  if ~isdir(dirName), dirName = fileparts(dirName); end
  if isempty(dirData), return; end
  dirIndex = [dirData.isdir];             % Find the index for directories
  fileList = {dirData(~dirIndex).name}';  % Get a list of the files
  if ~isempty(fileList)
    % exact search
    index = find(strcmp(File, fileList));
    if ~isempty(index)
      fileList = fileList(index);
      for i=1:numel(fileList)
        fileList{i} = fullfile(dirName,fileList{i});  % get the full path/file name
      end
      return
    end
    % more relaxed search
    index = find(~cellfun(@isempty,strfind(fileList, File)));
    if ~isempty(index)
      fileList = fileList(index);
      for i=1:numel(fileList)
        fileList{i} = fullfile(dirName,fileList{i});  % get the full path/file name
      end
      return
    else
      fileList = {};
    end
  end
  subDirs = {dirData(dirIndex).name};          % Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'}) ...
    & ~strncmp(subDirs, '.',1);                % Find index of subdirectories
                                               %   that are not '.' or '..' or hidden

  if recursive
    for iDir = find(validIndex)                  % Loop over valid subdirectories
      nextDir = fullfile(dirName,subDirs{iDir}); % Get the subdirectory path
      fileList = getAllFiles(nextDir, File);     % Recursively call getAllFiles
      if ~isempty(fileList), return; end
    end
  end

% ------------------------------------------------------------------------------
function [info, exe] = instrument_get_info(executable)
  % calls the instrument with --info or --help to attempt to get information
  % parse the output and return a structure.

  info = ''; exe = '';
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ; DISPLAY= ; ';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; DISPLAY= ; ';
  else           precmd = ''; end

  [p,f,e] = fileparts(executable);

  for name={ fullfile(p,f), f,  fullfile('.',f) }
    for ext={'','.out','.exe'}
      if isempty(dir([ name{1} ext{1} ])), continue; end
      for opt={' --info',' --help',' -h'}
        % look for executable and test with various extensions
        exe = [ name{1} ext{1} ];
        [status, result] = system([ precmd exe opt{1} ]);
        if any(status == [ -1 0 255 ])
          info = result;
          return
        end
      end
    end
  end

% ------------------------------------------------------------------------------
function result = instrument_compile(options)

  % compile the instrument using McCode
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ; DISPLAY= ; ';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; DISPLAY= ; ';
  else           precmd = ''; end

  if isempty(options.mccode)
      status=1;
      result = [ mfilename ': McCode (mcrun/mxrun) not found. Try "mccode check" first.' ];
  else
      disp([ mfilename ': Compiling instrument from ' options.instrument ...
        ' using ' options.mccode]);
      % assemble the command line: compile, no particle generated
      cmd = [ options.mccode ' --force-compile ' options.instrument ' --ncount=0' ];
      if isfield(options,'mpi') && options.mpi > 1
        cmd = [ cmd ' --mpi=1' ];
      end
      disp(cmd)
      p=pwd;
      [status, result] = system([ precmd cmd ]);
      cd(p)
  end

  % stop if compilation fails...
  if ~any(status == [ -1 0 255 ])
    disp(result);
    error([ mfilename ': ERROR: failed compilation of ' options.instrument ]);
  end

% ------------------------------------------------------------------------------
function [dynamic, static] = instrument_get_parameters(info)
  % analyze the --help or --info string and search for parameters with default values

  % we search
  % line  'Parameters: par1(type) par2(type) ...
  % lines 'Param: <name>=<val>  only for parameters which have default values
  % lines 'par (type)'
  % lines 'par (type) [default='val']'
  %
  % 'type' can be: double, string, int

  if isempty(info), return; end
  info = textscan(info, '%s','Delimiter','\n\r');
  info = info{1};

  dynamic = struct();
  static  = struct();

  % will set flag to true when we start to scan 'Instrument parameters'
  flag_instr_pars = false;

  % now we scan lines
  for l_index = 1:numel(info)
    this_line = strtrim(info{l_index});
    if strncmpi(this_line, 'Parameters:', length('Parameters:'))
      this_line = strtrim(this_line(length('Parameters: '):end));
      % regular expression to search for par(type)
      % store all these parameter names: double as variable, others as static
      % default value is set to []
      par=regexp(this_line,'(\w+)\s*\(double\)',    'tokens');
      for index=1:numel(par)
          name = par{index};
          if ~isempty(name), dynamic.(name{end}) = []; end
      end
      par=regexp(this_line,'(\w*)\s*\((string|int)\)','tokens');
      for index=1:numel(par)
          name = par{index};
          if ~isempty(name), static.(name{1}) = []; end
      end

    elseif strncmpi(this_line, 'Param:', length('Param:'))
      par = str2struct(this_line); c={}; f={};
      while isstruct(par)
        c = struct2cell(par); % get the content
        f = fieldnames(par);
        if isstruct(c{end}), par = c{end};
        else par = c; end
      end
      % assign the parameter depending on its type
      if ~isempty(c)
        if   isnumeric(c{end}), dynamic.(f{end}) = c{end};
        else                    static.(f{end})  = c{end}; end
      end
    elseif strncmpi(this_line, 'Instrument parameters', length('Instrument parameters'))
      flag_instr_pars = true;
    elseif flag_instr_pars
      % search for parameter with default value
      par=regexp(this_line,'(\w)+\s*\((double|string|int)\)\s*\[default=''(.*)''\]','tokens');
      if ~isempty(par)
          % tok has 3 elements: name, type, value
          for index=1:numel(par)
            var = par{index};
            name = var{1};
            val  = var{3};
            if ~isempty(str2num(val)), val = str2num(val); end
            if strcmp(var{2},'double'), dynamic.(name) = val;
            else                 static.(name)  = val; end
          end
          continue
      end
      % search for parameter without default value
      par=regexp(this_line,'(\w)+\s*\((double|string|int)\)','tokens');
      % tok has 2 elements: name, type
      for index=1:numel(par)
        var = par{index};
        name = var{1};
        if strcmp(var{2},'double'), dynamic.(name) = [];
        else                 static.(name)  = []; end
      end
    end
  end