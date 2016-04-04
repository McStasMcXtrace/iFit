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

% this function is called when the MeX is not present/compiled -> use Binary
% the result is not interpreted, i.e. the actual data is usually written to
% an external file that should be read afterwards.

result = [];

% use looktxt bin when available -----------------------------------------------
if isunix, precmd = 'LD_LIBRARY_PATH= ; '; else precmd=''; end

% handle input arguments

% assemble command line for the binary call
this_path = fileparts(which(mfilename));
cmd       = fullfile(this_path, mfilename);
if ispc, cmd=[ cmd '.exe' ]; end
cmd       = [ cmd ' ' sprintf('%s ', varargin{:}) ];
disp(cmd)

% launch the command
[status, result] = system([ precmd cmd ]);


