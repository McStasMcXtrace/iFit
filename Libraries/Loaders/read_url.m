function filename = read_url(filename, option)
% READ_URL Reads a distant file and returns a local file for further processing.
%   Temporary files are NOT removed afterwards. Think about removing them
%   after further processing.
%
%   data=READ_URL(url) get a distant resource.
%
%   data=READ_URL(url, tmpdir) extract compressed data from given file, 
%   in given temporary directory. It is recommended to use a RAMdisk for higher
%   efficiency, such as '/dev/shm' or '/run/shm'.
%
% Input:  filename: URL starting with file://, ftp:// http:// and https://
% output: local file path
% Example: y=read_url('http://www.jcamp-dx.org/lisms/cc/bruker1.dx'); delete(y); ischar(y)
%
% 
% See also: read_compressed, iLoad

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    s.name     ='URL';
    s.method   =mfilename;
    s.extension={'url'};
    return
end

if nargin < 2, option = ''; end

if isempty(option), option= tempdir; end
if ischar(option) && ~isempty(dir(option)) && ~isdir(option)
  option = fullparts(option); % get path
end

% handle file on the internet
if strncmp(filename, 'http://', length('http://')) ...
 | strncmp(filename, 'https://',length('https://')) ...
 | strncmp(filename, 'ftp://',  length('ftp://'))
  tmpfile = tempname(option);
  % Keep file extension, may be useful for iData load
  [filepath,name,ext] = fileparts(filename);
  tmpfile = [tmpfile ext];
  use_wget = false;
  if ~usejava('jvm')
    use_wget = true;
  else
    % access the net. Proxy settings must be set (if any).
    try
      % write to temporary file
      tmpfile = urlwrite(filename, tmpfile);
    catch ME
      if config.verbosity
        disp(ME.message);
      end
      use_wget = true;
    end
  end
  if use_wget
    % Fall back to using wget
    cmd = ['wget ' filename ' -O ' tmpfile]; disp(cmd)
    [status, result] = system(cmd);
    if status
      if config.verbosity
        disp(result);
      end
      error([ mfilename ': Can not get URL ' filename ]);
    end
  end
  filename = tmpfile;
end

if strncmp(filename, 'file://', length('file://'))
  filename = filename(7:end); % remove 'file://' from local name
end
