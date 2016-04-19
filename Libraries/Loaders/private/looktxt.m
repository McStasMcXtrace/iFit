function [result,status] = looktxt(varargin)
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
persistent target

result = []; status = 127;

% use looktxt bin when available -----------------------------------------------

% handle input arguments

% assemble command line for the binary call

% identify if we use a local or global (system) looktxt 

% local looktxt
this_path = fileparts(which(mfilename));
cmd       = fullfile(this_path, mfilename);
if ispc, ext = '.exe'; else ext = ''; end

if isempty(target)
  % we test if the executable files exist
  % test in order: global(system), local, local_arch
  for try_target = { [ cmd '_' computer('arch') ext ], [ cmd ext ], [ 'looktxt' ext ]}
      [status, result] = system(try_target{1});
      if status == 0
          target = try_target{1};
          disp([ mfilename ': Bin is valid from ' target ]);
          break
      end
  end
  if isempty(target)
      error([ mfilename ': no valid executable file found.' ])
      return
  end
end
cmd       = [ target ' ' sprintf('%s ', varargin{:}) ];
disp(cmd)

% launch the command
[status, result] = system([ cmd ]);


