function s = read_anytext(varargin)
% import any text using 'looktxt'.
%
% the importation consists in performing the folling tasks:
% * handle arguments, looking for options (possibly with "string") and filenames
% * launch looktxt as MeX and MATfile format on temporary file
% * import the MAT file as a structure

% handle input arguments =======================================================
if isempty(varargin)
  looktxt('--help');
  return
end

argv={};
format = 'MATfile'; % default output format, but can be overriden
user.outfile = '';
user.format  = '';
remove_tempname = 0;

% split the arguments, but retain those that contain "" and '' quotes
for index=1:length(varargin)
  arg = varargin{index};
  if ~ischar(arg)
    fprintf(1, 'read_anytext: argument is of class %s. Only char allowed. Ignoring\n', ...
      class(arg));
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

% launch looktxt with MATfile format and temporary file name ===================

% specify default output file name and format (if not defined by user)
if isnumeric(user.outfile) && user.outfile <= length(argv)
  user.outfile = argv{user.outfile}; 
end
if isnumeric(user.format) && user.format <= length(argv)
  user.format = argv{user.format}; 
end

if isempty(user.outfile)
  user.outfile     = [ tempname '.mat' ];  % usually in TMP directory
  argv{end+1} = [ '--outfile=' user.outfile ];
  remove_tempname = 1;
elseif strncmp(user.outfile, '--outfile=', length('--outfile='))
  user.outfile= user.outfile((length('--outfile=')+1):end);
end
if isempty(user.format)
  user.format =          'MATfile';
  argv{end+1} = '--format=MATfile';
elseif strncmp(user.format, '--format=', length('--format='))
  user.format= user.format((length('--format=')+1):end);
end

% call looktxt >>>>
s = looktxt(argv{:});

% import the MAT file from the temporary file, into structure ==================
if strcmp(format, 'MATfile')
  s = user.outfile;
  try
    s = load(user.outfile);
    % check if there is only one struct field at first level then access it (probably
    % the temporary variable name).
    f = fieldnames(s);
    if length(f) == 1
      s = s.(f{1});
    end
  catch
    s= user.outfile;
  end
end

% delete temporary file ========================================================
if remove_tempname
  delete(user.outfile);
end
