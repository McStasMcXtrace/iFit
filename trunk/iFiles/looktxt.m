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
% Example: looktxt -c -s PARAM -s DATA filename
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
% looktxt  version 1.0.0 (2 Sept 2006) by Farhi E. [farhi@ill.fr]

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
  looktxt_exe = [ fileparts(which('looktxt')) filesep looktxt_exe ];
  exec = [ looktxt_exe ' --outfile=' file ' ' args ];
  [s,w] = system(exec);
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
else
  warning('iLoad:looktxt','looktxt: can not find executable. Re-install/compile looktxt');
end
cd(p);
