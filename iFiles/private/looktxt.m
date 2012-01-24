% [structs] = looktxt('filename', '[options]') Import text data
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
%
% Usual options are: --fast --fortran --binary --force --catenate --comment=NULL
% List of all options can be obtained using: looktxt --help
%
% looktxt  version 1.1 $Revision: 1.2 $ (24 Sept 2009) by Farhi E. [farhi@ill.fr]

% if we come there, that's because the mex file is not compiled.
% we first try to install it, and if it fails, we go for the CC version

% Looktxt: Search and export numerics in a text/ascii file
% Copyright (C) 2009  E. Farhi, Institut Laue Langevin
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%


function data = looktxt(args)
data = [];
if nargin == 0, args=''; end

if exist('mex') && exist('texmex.c')
  looktxtmex = which('texmex.c');
  looktxtpath = fileparts(looktxtmex);
  exec = [ 'mex -O -output ' looktxtpath filesep 'looktxt ' looktxtmex ];
  disp(exec);
  try
    eval(exec);
  catch
    if ispc
      % assume we use LCC
      exec=['mex -O -v -output ' looktxtpath filesep 'looktxt ' looktxtmex ' -L"' fullfile(matlabroot,'sys','lcc','lib') '" -lcrtdll' ];
      disp(exec);
      eval(exec);
    else
      error([ mfilename ': Installation failed. Check your C compiler installation.' ]);
    end
  end
  rehash
  if nargin > 0, data = looktxt(args); end
  return
end

% if we come here, mex fails, and we try the CC version
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

% automatic compilation of looktxt using CC
if ~exist([ file '.m' ],'file') & s ~= 0  % executable not found ?
    
  % first check if system wide looktxt exists
  ok=0;
  [s,w1] = system([ looktxt_exe ]);
  if (s==0) ok=1; end  % looktxt is installed.
  
  % then check if local looktxt exists
  if ~ok
    looktxtc=which('looktxt.c');  % where C code is
    pathstr = fileparts(looktxtc);
    [s,w2] = system(fullfile(pathstr, looktxt_exe));
    if (s==0) ok=1; end  % looktxt is installed (local)
  end
  
  if ok
    disp(w)
    return
  end
  
  
  % then tries to install looktxt
  if isempty(looktxtc), error('Can not install looktxt as source code is unavailable'); end
  
  % identify C compiler
  cc     = getenv('CC');     
  if isempty(cc),
    [s,w] = system('gcc'); 
    if s == 0, cc = 'gcc';
    else cc='cc'; end
  end
  cflags = getenv('CFLAGS'); if isempty(cflags), cflags = '-O2'; end
  
  % build executable
  warning('looktxt: can not find executable. Attempting to re-install/compile looktxt');
  disp([ cc ' ' cflags ' -o ' pathstr filesep looktxt_exe ' ' looktxtc ]);
  [s,w] = system([ cc ' ' cflags ' -o ' pathstr filesep looktxt_exe ' ' looktxtc ]);
  if s ~= 0
    disp(w);
    return
  end
  
  disp('looktxt: Testing the validity of executable')
  [s,w] = system(fullfile(pathstr, looktxt_exe));
  if s ~= 0
    disp(w);
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
