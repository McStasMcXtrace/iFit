% [structs] = looktxt('filename [options]') Import text data
% Usage: looktxt [options] file1 file2 ...
% Action: Search and export numerics in a text/ascii file.
%    This program analyses files looking for numeric parts
%    Each identified numeric field is named and exported
%    into a structure with fields
%      ROOT.SECTION.FIELD = VALUE
%    In order to sort your data, you may specify as many --section
%    and --metadata options as necessary
% Return: a single structure or a cell of structures
%    If the structure can not be evaluated, the raw Matlab script is returned.
% Example:   looktxt -c -s PARAM -s DATA filename
% Usual options are: --fast --fortran --binary --force --catenate --comment=NULL
%
% Useful options when used from Matlab:
% --binary   or -b    Stores numerical matrices into an additional binary
%                     float file, which makes further import much faster.
% --catenate or -c    Catenates similar numerical fields
% --force    or -F    Overwrites existing files
% --fortran --wrapped Catenates single Fortran-style output lines with
%                     previous matrices
% --headers  or -H    Extracts headers for each numerical field
% --section=SEC       Classifies fields into section matching word SEC
%       -s SEC
% --metadata=META     Extracts lines containing word META as user meta data
%         -m META
% --fast              Uses a faster reading method, requiring numerics
%                     to be separated by \n\t\r\f\v and spaces only
% --makerows=NAME     All fields matching NAME are transformed into row vectors
% List of all options can be obtained using: looktxt --help
%
% looktxt  version 1.0.7 (10 July 2009) by Farhi E. [farhi@ill.fr]

function data = looktxt(args)
data = [];
if nargin == 0, args=''; end
if iscellstr(args)
  data = cell(length(args), 1);
  for index=1:length(args)
    data{index} = looktxt(args{index});
  end
  return
end
% handle single file name
file = tempname;  % results will go there
if ispc, looktxt_exe = 'looktxt.exe';
else     looktxt_exe = 'looktxt';
end
% execute command
exec = [ looktxt_exe ' --outfile=' file ' ' args ];
disp(exec);
[s,w] = system(exec);
  
% check if result has been generated, else try again with local executable
if ~exist([ file '.m' ],'file')
  exec = [ fileparts(which('looktxt')) filesep looktxt_exe ' --outfile=' file ' ' args ];
  [s,w] = system(exec);
end
if ~exist([ file '.m' ],'file') & s ~= 0  % executable not found
  if isempty(args) || strfind(args, '--help') || strfind(args, '--version') 
    help looktxt
    return
  end
  warning('looktxt: can not find executable. Attempting to re-install/compile looktxt');
  cc     = getenv('CC');     
  if isempty(cc),
    [s,w] = system('gcc'); 
    if s == 0, cc = 'gcc';
    else cc='cc'; end
  end
  cflags = getenv('CFLAGS'); if isempty(cflags), cflags = '-O2'; end
  looktxtc=which('looktxt.c');  % where C code is
  if isempty(looktxtc), error('Can not install looktxt as source code is unavailable'); end
  path = fileparts(looktxtc);
  disp([ cc ' ' cflags ' -o ' path filesep looktxt_exe ' ' looktxtc ]);
  [s,w] = system([ cc ' ' cflags ' -o ' path filesep looktxt_exe ' ' looktxtc ]);
  disp('looktxt: Testing the validity of executable')
  [s,w] = system([ path filesep looktxt_exe ]);
  if s ~= 0
    error('looktxt: Failed to install/compile looktxt. Please install it manually.');
  else
    disp('looktxt: OK, executable is functional');
    % now re-try with executable
    [s,w] = system(exec);
  end
end
disp(w);
% now import structure
p=pwd;
if exist([ file '.m' ],'file')
  cd(tempdir);
  try
    [pathstr, name, ext] = fileparts(file);
    data = feval(name);
  catch
    disp('looktxt: error evaluating result');
    lasterr
  end
  try
    % remove tmp filename
    delete([ file '.m' ]);
    if ~isempty(strfind(args, '--binary')) | ~isempty(strfind(args, '-b'))
      delete([ file '.bin' ]);
    end
  catch
  end
elseif (s ~= 0)
    error([ 'looktxt: Failed to import ' args ]);
end
cd(p);
