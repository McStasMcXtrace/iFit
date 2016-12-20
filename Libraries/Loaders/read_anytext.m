function s = read_anytext(varargin)
% import any text using 'looktxt'.
%   data = read_anytext(filename, options...)
%
% The possible options are:
%'--catenate'	 Catenates similar numerical fields (which have similar dimensions 
%                and names. Recommended.
%'--fast'      When numerical data blocks only use isspace(3) separators 
%                (\n \r \f \t \v and space), the reading can be made faster with 
%                even lower memory requirements. Recommended.
%'--headers'   Extracts headers for each numerical field. Recommended.
%'--wrapped'   Catenates single wrapped output lines with previous matrices 
%                (e.g. caused by the 80 chars per line limit in old data formats 
%                written by fortran codes). Recommended.
%'--section=SEC' Classifies fields into sections matching word SEC. This option 
%                can be repeated with different SEC words.
%'--metadata=META' Extracts lines containing word META as user metadata. This 
%                option can be repeated with different META items.
%'--makerows=NAME'	When a numerical data block label matching NAME is found, it 
%                is transformed into a row vector. This may be used for wrapped 
%                files (--wrapped option). This option can be repeatedwith as 
%                different NAME tokens.
%'--help'      Lists all possible options.
%'--silent'    Suppress processing messages except errors.
%
% the importation consists in performing the following tasks:
% * handle arguments, looking for options (possibly with "string") and filenames
% * launch looktxt with Matlab/binary format on temporary file
% * import the MAT file as a structure
%
% read_anytext('compile') check looktxt and recompiles if needed
% (c) E.Farhi, ILL. License: EUPL.

% we choose NOT to use the looktxt mex file due to SEGV under Matlab.

persistent compiled

s = [];

% the configuration must now be a string: derive executable and output format
% the default is  Matlab binary, which is fast and well supported in all cases
output     = 'Matlab'; 

% *** test executable ==========================================================
if ~isempty(varargin) && strcmp(varargin{1}, 'compile'), compiled = ''; end

% test if binary is requested and exists.
if isempty(compiled) % only the first time it starts or when explicitly requested
  
  % test if bin is requested and exists, else compiles
  if ~isempty(varargin) && strcmp(varargin{1}, 'compile')
       compiled = read_anytext_compile_binary('compile'); % force
  else compiled = read_anytext_compile_binary; end
  if ~isempty(compiled), 
    config='bin'; 
  end
  if isempty(compiled)
    error('%s: ERROR: Can''t compile looktxt executable (Binary).', ...
          mfilename);
  end
  
  if ~isempty(varargin) && strcmp(varargin{1}, 'compile')
    varargin(1) = [];
  end
end

% *** handle input arguments ===================================================

if isempty(varargin)
  s = looktxt('--help');
  return
elseif strcmp(varargin{1}, 'compile') || strcmp(varargin{1}, 'check')
  s = compiled;
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
if isempty(user.format) && ~isempty(output), user.format = output; end

% when MATFile or Matlab and no output file set, use temporary
if isempty(user.outfile)
  if strcmp(user.format, 'MATFile')
    user.outfile= [ tempname '.mat' ];  % usually in TMP directory
  elseif strcmp(user.format, 'Matlab')
    user.outfile= [ tempname '.m' ];    % usually in TMP directory
  end
  remove_tempname = 1;
end

if strcmp(user.format, 'Matlab')
  argv{end+1} = '--names_root=data';  % for easier evaluation in read_anytext_eval
  argv{end+1} = '--binary';           % faster
end

% send options as looktxt arguments
if ~isempty(user.format)
  argv{end+1} = [ '--format=' user.format ];
end
if ~isempty(user.outfile)
  argv{end+1} = [ '--outfile=' user.outfile ];
end

% the filename should be the first argument. Test if binary
[fid,message] = fopen(argv{1}, 'r');
if fid ~= -1
  file_start = fread(fid, 1000, 'uint8=>char')';
  fclose(fid);
  if length(find(file_start >= 32 & file_start < 127))/length(file_start) < 0.4,
    return  % this is a binary file. Skip.
  end
end

% *** call looktxt =============================================================
s = [];
[result,status] = looktxt(argv{:}); % send to looktxt.m to launch bin

% *** import the data (user.format) ============================================
if strcmp(user.format, 'MATFile')
  % import the MAT file from the temporary file, into structure 
  if ~isempty(user.outfile) && ischar(user.outfile) && ~isempty(dir(user.outfile))
    s = load(user.outfile); % must be a MAT-file
  end
elseif strcmp(user.format, 'Matlab')
   % import the Matlab script file from the temporary file, into structure 
  if ~isempty(user.outfile) && ischar(user.outfile) && ~isempty(dir(user.outfile))
    s = read_anytext_eval(user.outfile);
  end
end
if isstruct(s)
  % check if there is only one struct field at first level then access it (probably
  % the temporary variable name).
  f = fieldnames(s);
  if length(f) == 1
    s = s.(f{1});
  end
else % generated data is empty ??!! 
  fprintf(1, [ '%s: ERROR: empty generated data %s: looktxt %s\n' ], ...
  mfilename, user.format, sprintf('%s ', argv{:}));
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
    s.Data.read_anytext = { compiled; user.format; [ 'looktxt ' sprintf('%s ', argv{:}) ] };
  end
end

% -------------------------------------------------------------------------       
function d=read_anytext_eval(str)
  % evaluate a Matlab/bin generated output in a reduced environment
  d = [];
  if ~isdeployed
    run(str);  % root level is 'data' as set with root_name=data
  else
    % must remove 'function' lines and use the 'bin_ref' function below
    % as 'function' is not allowed in standalone.
    s = fileread(str);
    if ~isempty(s)
      % search '%% To import' which is generated at start of file
      pos1 = strfind(s, '%% To import'); if isempty(pos1), pos1=1; end
      % search '%% in-line function to read binary blocks' generated at end
      pos2 = strfind(s, '%% in-line function to read binary blocks');
      if isempty(pos2), pos2=numel(s); end
      f    = s(pos1:pos2); 
      eval(f); % evaluate cleaned file, creates 'data' structure
    end
  end
  % get the result of evaluation
  if     exist('data','var'), d=data; 
  elseif exist('ans','var'),  d=ans; end
  
function d=bin_ref(f,b,m,n)
  [fid,mess]=fopen(f,'rb');
  if fid == -1, disp([ 'Error opening bin file ' f ': ' mess ]); d=[]; return; end
  fseek(fid,b,-1);
  d=fread(fid,m*n,'double'); fclose(fid);
  if m*n ~= numel(d), 
    disp([ 'File ' f ': read ' num2str(numel(d)) ' elements but expected ' mat2str([ m n ]) ]); 
    f=dir(f); disp(f); 
  end
  d=reshape(d,n,m);
  d=d'; return




% -------------------------------------------------------------------------
function compiled = read_anytext_compile_binary(compile)
  % compile looktxt as binary when does not exist yet
  
  compiled = ''; 
  if nargin == 0, compile = ''; end
  if ismac,  precmd = 'DYLD_LIBRARY_PATH= ;';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
  else precmd=''; end
  
  if ispc, ext='.exe'; else ext=''; end
  this_path = fullfile(fileparts(which(mfilename)),'private');
  
  % test if we still have a looktxt MeX: remove it as it will probably give SEGV
  % check if it exists and is valid
  if exist(which('looktxt')) == 3
    disp([ mfilename ': WARNING: Removing conflicting MeX file ' which('looktxt') ])
    delete(which('looktxt')); compiled = '';
  end
  
  % try in order: global(system), local, local_arch
  for try_target={ ...
          fullfile(this_path, [ 'looktxt_' computer('arch') ext ]), ...
          fullfile(this_path, [ 'looktxt' ext ]), ...
          [ 'looktxt' ext ], 'looktxt'}
      
    [status, result] = system([ try_target{1} ' --help' ]); % run from Matlab

    if status == 0 && nargin == 0
        % the executable is already there. No need to make it .
        compiled = try_target{1};
        return
    end
  end
  
  % when we get there, compile looktxt_arch, not existing yet
  target = fullfile(this_path, [ 'looktxt_' computer('arch') ext ]);

  % search for a C compiler
  cc = '';
  for try_cc={getenv('CC'),'cc','gcc','ifc','pgcc','clang','tcc'}
    if ~isempty(try_cc{1})
      [status, result] = system([ precmd try_cc{1} ]);
      if status == 4 || ~isempty(strfind(result,'no input file'))
        cc = try_cc{1};
        break;
      end
    end
  end
  if isempty(cc)
    if ~ispc
      disp([ mfilename ': ERROR: C compiler is not available from PATH:' ])
      disp(getenv('PATH'))
      disp([ mfilename ': You may have to extend the PATH with e.g.' ])
      disp('setenv(''PATH'', [getenv(''PATH'') '':/usr/local/bin'' '':/usr/bin'' '':/usr/share/bin'' ]);');
    end
    error('%s: Can''t find a valid C compiler. Install any of: gcc, ifc, pgcc, clang, tcc\n', ...
    mfilename);
  else
    try
      fprintf(1, '%s: compiling looktxt binary (using %s)...\n', mfilename, cc);
      cmd={cc, '-O2','-o',target, ...
        fullfile(this_path,'looktxt.c'),'-lm'};
      cmd = sprintf('%s ',cmd{:});
      disp(cmd)
      [status, result] = system([ precmd cmd ]);
      if status == 0
        compiled = target;
      end
    end
  end

  if isempty(compiled) && ~isempty(compile)
    error('%s: Can''t compile looktxt.c binary\n       in %s\n', ...
        mfilename, fullfile(this_path));
  end
  

