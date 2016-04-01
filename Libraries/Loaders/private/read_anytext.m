function s = read_anytext(varargin)
% import any text using 'looktxt'.
%
% the importation consists in performing the following tasks:
% * handle arguments, looking for options (possibly with "string") and filenames
% * launch looktxt as MeX and MATfile format on temporary file
% * import the MAT file as a structure
%
% read_anytext('compile') creates looktxt MeX or binary
% read_anytext('config')  update  iLoad config

% different modes for execution:
%  The configuration is extracted from iLoad('config')
% Executable      Output
%   MeX             Mem               direct memory allocation
%   Mex             Matlab (m+bin)    a Matlab script+binary file
%   MeX             MATFile           a MAT file
%   Bin             Matlab (m+bin)
%   Bin             MATFile

persistent config
persistent compiled

if nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1},'config')
  config = [];
end

% get configuration to use for 'looktxt'
if isempty(config)
  if exist('iLoad')
    config  = iLoad('config');
    if nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1},'config')
      return
    end
  else
    config.MeX = 'default';
  end
end
if isfield(config, 'MeX')
  tmp = config.MeX; config=[]; % to avoid warning when erasing structure
  config = tmp;
end

% *** determine executable (MeX/binary) and output format (mem/mat/matlab) =====
if isstruct(config) 
  if isfield(config, 'looktxt')
    tmp = config.looktxt;
  elseif isfield(config, 'read_anytext')
    tmp = config.read_anytext;
  end
  config = tmp;
end

if isempty(config) || ~ischar(config)
  disp([ mfilename ': WARNING: invalid text importer configuration. Using ''default''.' ]); 
  disp(config);
  config = 'default'; 
end

% the configuration must now be a string: derive executable and output format
config = lower(config);
executable = '';
output     = '';

% we interpret the requested configuration
if ~isempty(strfind(config, 'mex')) || ~isempty(strfind(config, 'mem')), executable = 'MeX'; output='MeX'; end
if ~isempty(strfind(config, 'bin')),     executable = 'Binary'; output='Matlab'; end
if ~isempty(strfind(config, 'matlab')),  output='Matlab'; 
elseif ~isempty(strfind(config, 'mat')), output='MATFile'; end

if ~isempty(strfind(config, 'default')) || ~isempty(strfind(config, 'auto')) ...
  || isempty(executable) || isempty(output) % default/auto choices
  % default mex/bin choice set by the system type: same as config.MeX='default'
  if isempty(executable), executable = 'mex'; end
  if isempty(output)
    if ispc || ismac, output = 'MeX';         % in memory
    else              output = 'MATFile'; end % Linux: avoid MeX/mem which may be unstable (SEGV)
  end  
end

% *** test executable ==========================================================
if ~isempty(varargin) && strcmp(varargin{1}, 'compile'), compiled = 0; end
executable = lower(executable);

% test if mex is requested and exists. When fails -> try bin.
if isempty(compiled) || ~compiled % only the first time it starts or when explicitly requested
  if ~isempty(strfind(executable, 'mex'))
    if ~isempty(varargin) && strcmp(varargin{1}, 'compile')
      s = read_anytext_compile_mex('compile'); % force
    else s = read_anytext_compile_mex; end
    if isempty(s), executable = 'bin'; % will try bin
    end
  end
  
  % test if bin is requested and exists, else compiles
  if ~isempty(strfind(executable, 'bin'))
    if ~isempty(varargin) && strcmp(varargin{1}, 'compile')
      s = read_anytext_compile_binary('compile'); % force
    else s = read_anytext_compile_binary; end
    if isempty(s), 
      error('%s: ERROR: Can''t compile looktxt executable (MeX nor Binary).', ...
          mfilename); 
    end
  end
  compiled=1;
  
  if ~isempty(varargin) && strcmp(varargin{1}, 'compile')
    varargin(1) = [];
  end
end

% *** handle input arguments ===================================================
s = [];
if isempty(varargin)
  s = looktxt('--help');
  return
end

argv={};

user.outfile = '';
user.format  = '';
remove_tempname = 0;

% split the arguments, but retain those that contain "" and '' quotes
for index=1:length(varargin)
  arg = varargin{index};
  if ~ischar(arg)
    fprintf(1, '%s: WARNING: argument %i is of class %s. Only char allowed. Ignoring\n', ...
      mfilename, index, class(arg));
  else
    split = strread(arg,'%s','delimiter',' ;'); % split argument with common delimiters
    i_split = 1;
    while i_split <= length(split)
      this_split = split{i_split}; c='';
      if any(this_split == '''')  % is there an argument which contains a string ?
        c = '''';
      elseif any(this_split == '"')
        c = '"';
      end
      while ~isempty(c) && i_split < length(split)
        % assemble arguments until the last one that ends with the same delimiter
        i_split = i_split+1;
        this_split = [ this_split ' ' split{i_split} ];
        if any( split{i_split} == c )
          c = ''; % we have found the second string delimiter: end while loop
        end
      end
      if strncmp(this_split, '--outfile=', length('--outfile=')) ...
        user.outfile = length(argv)+1;
      elseif any(strcmp( this_split, {'-o','--outfile'})) % file name follows -o
        user.outfile = length(argv)+2;
      end
      if strncmp(this_split, '--format=', length('--format=')) ...
        user.format  = length(argv)+1;
      elseif any(strcmp( this_split, {'-f','--format'})) % format follows -f 
        user.format  = length(argv)+2;
      end
      argv{end+1} = this_split; % store this bit as an argument for looktxt
      i_split = i_split+1;
    end % for
  end
end

% clean-up format and outfile options
if isnumeric(user.outfile) && user.outfile <= length(argv)
  user.outfile = argv{user.outfile}; 
end
if isnumeric(user.format) && user.format <= length(argv)
  user.format = argv{user.format};
end
if strncmp(user.outfile, '--outfile=', length('--outfile='))
  user.outfile= user.outfile((length('--outfile=')+1):end);
end
if strncmp(user.format, '--format=', length('--format='))
  user.format = user.format((length('--format=')+1):end);
end

% when no format specified use config
if isempty(user.format)
  switch lower(output)
  case 'mex'
    user.format = 'MeX';
  case {'mat','matfile'}
    if exist('looktxt') == 3 % we make sure the MATFile is used with the MeX
      user.format = 'MATFile';
    else
      user.format = 'Matlab';
    end
  otherwise
    user.format = 'Matlab';
  end
end

% when MATFile or Matlab and no output file set, use temporary
if isempty(user.outfile)
  if strcmp(user.format, 'MATFile')
    user.outfile= [ tempname '.mat' ];  % usually in TMP directory
  elseif strcmp(user.format, 'Matlab')
    user.outfile= [ tempname '.m' ];  % usually in TMP directory
    argv{end+1} = '--binary';
  end
  remove_tempname = 1;
end

% send options as looktxt arguments
if ~isempty(user.format)
  argv{end+1} = [ '--format=' user.format ];
end
if ~isempty(user.outfile)
  argv{end+1} = [ '--outfile=' user.outfile ];
end

% the filename should be the first argument
[fid,message] = fopen(argv{1}, 'r');
if fid ~= -1
  file_start = fread(fid, 1000, 'uint8=>char')';
  if length(find(file_start >= 32 & file_start < 127))/length(file_start) < 0.4,
    return  % this is a binary file. Skip.
  end
end

% *** call looktxt =============================================================
if strcmp(executable, 'mex')
  % pure MEX call. No temporary file. May cause SEGV. faster by 15%.
  s = looktxt(argv{:});
  result = ''; status=0;
elseif strcmp(executable, 'bin')
  s = [];
  [status,result] = looktxt(argv{:}); % send to looktxt.m to launch bin
end

% *** import the data (user.format) ============================================
if strcmp(user.format, 'MATFile')
  % import the MAT file from the temporary file, into structure 
  if ~isempty(user.outfile) && ischar(user.outfile) && ~isempty(dir(user.outfile))
    try
      s = load(user.outfile); % must be a MAT-file
      % check if there is only one struct field at first level then access it (probably
      % the temporary variable name).
      f = fieldnames(s);
      if length(f) == 1
        s = s.(f{1});
      end
    catch
      fprintf(1, [ '%s: ERROR: looktxt ' argv{:} '\n' ], mfilename);
    end
  end
elseif strcmp(user.format, 'Matlab')
   % import the Matlab script file from the temporary file, into structure 
  if ~isempty(user.outfile) && ischar(user.outfile) && ~isempty(dir(user.outfile))
    try
      run(user.outfile);
      s=ans;
      % check if there is only one struct field at first level then access it (probably
      % the temporary variable name).
      f = fieldnames(s);
      if length(f) == 1
        s = s.(f{1});
      end
    end
  end
end

% delete temporary file
if remove_tempname && ~isempty(dir(user.outfile))
  % delete any [user.outfile].* file
  [p,f] = fileparts(user.outfile);
  delete(fullfile(p,[ f '.*' ]));
end

% convert the Headers field into Attributes
if isstruct(s)
  if isfield(s, 'Headers')
    s.Attributes = s.Headers;
    s=rmfield(s, 'Headers');
  end

  s=orderfields(s);

  if isfield(s, 'Data') && isstruct(s.Data)
    s.Data = orderfields(s.Data);
  end
  if isstruct(s.Data) && ~isfield(s.Data,'read_anytext')
    s.Data.read_anytext = { executable, output, [ 'looktxt ' sprintf('%s ', argv{:}) ] };
  end
end

% ------------------------------------------------------------------------------
function compiled = read_anytext_compile_mex(compile)
  % compile looktxt as MeX
  
  compiled = '';
  if isdeployed, return; end
  
  % check if it exists and is valid
  if exist(which('looktxt')) == 3
    try
      compiled = looktxt('--version');
      compiled = which('looktxt');
    end
  end
  if ~isempty(compiled) && nargin == 0, return; end
  
  this_path = fileparts(which(mfilename));
  % attempt to compile MeX
  fprintf(1, '%s: compiling looktxt mex...\n', mfilename);
  try
    cmd={'-O','-output',fullfile(this_path,'looktxt'), ...
         fullfile(this_path,'looktxt.c'),'-DUSE_MEX'};
    disp([ 'mex ' sprintf('%s ', cmd{:}) ]);
    mex(cmd{:});
    compiled = 'mex';
  catch
    try
      % Windows/LCC compiler support
      cmd={'-O','-output',fullfile(this_path,mfilename),...
           fullfile(this_path,'looktxt.c'),'-DUSE_MEX', ...
           ['-L"' fullfile(matlabroot,'sys','lcc','lib') '"' ],'-lcrtdll'};
      disp([ 'mex ' sprintf('%s ', cmd{:}) ]);
      mex (cmd{:});
      compiled = 'mex';
    catch
      error('%s: Can''t compile looktxt.c mex\n       in %s\n', ...
        mfilename, fullfile(this_path));
    end
  end

function compiled = read_anytext_compile_binary(compile)
  % compile looktxt as binary when does not exist yet
  
  compiled = '';
  
  if isdeployed, return; end
  if ispc, ext='.exe'; else ext=''; end
  this_path = fileparts(which(mfilename));
  % attempt to compile as binary, when it does not exist yet
  this_path = fileparts(which(mfilename));
  target    = fullfile(this_path, 'looktxt', ext);
  % launch the command
  [status, result] = system(target);
  if status == 0 && nargin == 0
    % the executable is already there. No need to make it .
    compiled = target; 
    return
  end

  fprintf(1, '%s: compiling looktxt binary...\n', mfilename);
  try
    cmd={'-f', fullfile(matlabroot,'bin','matopts.sh'), '-DUSE_MAT', ...
         '-O', '-output', target, ...
         fullfile(this_path,'looktxt.c'), '-lmat', '-lmx'};
    disp([ 'mex ' sprintf('%s ', cmd{:}) ]);
    mex(cmd{:});
    compiled = target;
  catch
    error('%s: Can''t compile looktxt.c binary\n       in %s\n', ...
        mfilename, fullfile(this_path));
  end

