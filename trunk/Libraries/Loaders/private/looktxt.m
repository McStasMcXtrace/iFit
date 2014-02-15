function result = looktxt(varargin)
% result = looktxt(file, options, ...)
%
% Action: Search and export numerics in a text/ascii file.
% This program analyses files looking for numeric parts
% Each identified numeric field is named and exported
% into an output filename.

%     export MATLABROOT=/opt/MATLAB/R2010a
%     export ARCH=glnxa64
%     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MATLABROOT/bin/$ARCH:$MATLABROOT/sys/os/$ARCH
%     gcc -I$MATLABROOT/extern/include -L$MATLABROOT/bin/$ARCH -DUSE_MAT -O2 -o looktxt -lmat -lmx looktxt.c

% this function is called when the MeX is not present/compiled

% first time, if the mex is not present, we try to compile it

persistent compiled

result = [];

% check if we use the looktxt executable, or need to compile the MeX (only once)
if ~isdeployed && (isempty(compiled) || strcmp(varargin{1}, 'compile'))
  compiled = read_anytext('compile');
  
  % succeed mex: rethrow looktxt call to the mex
  if strcmp(compiled, 'mex')
    rehash
    % rethrow command with new MeX
    clear functions
    if exist(mfilename) == 3
      result = looktxt(varargin{:});
    elseif nargin > 0
      error('%s: MeX compiled successfully. Rethrow command to use it.');
    end
  elseif ~strcmp(compiled, 'bin')
    % no bin, no mex : fail to compile and return
    return
  end
end

% use looktxt bin when available -----------------------------------------------
if ~strcmp(compiled, 'bin'), return; end
% handle input arguments

% assemble command line for the binary call
this_path = fileparts(which(mfilename));
cmd       = fullfile(this_path, mfilename);
cmd       = [ cmd ' ' sprintf('%s ', varargin{:}) ];
disp(cmd)

% launch the command
[status, result] = system(cmd);

% record if looktxt binary call works for the next time
if status ~= 0, compiled = 'no'; end


