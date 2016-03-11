function result = cif2hkl(varargin)
% result = cif2hkl(file_in, file_out, lambda, mode, verbose)
%
% Read a CIF/CFL/SHX/PCR crystallographic description
%       and generates a HKL F^2 reflection list.
%
% Input:
%   file_in:   name of CIF/CFL/PCR/ShellX file (char)
%   file_out:  name of output file, or empty for '<file_in>.laz' (char)
%   lambda:    minimum wavelength used for generation of HKL peaks. Default is 0.5  (scalar, Angs)
%   mode:      file conversion mode 'p'=powder, 'x'=single crystal, '-'=no/unactivated
%   verbose:   verbosity level, 0 or 1 (scalar)
% Output:
%   result:    command result (char)
%   file_out is written when mode is not '-'
%
% Revision: $date$

% this function is called when the MeX is not present/compiled

persistent compiled

result = [];
if isunix, precmd = 'LD_LIBRARY_PATH= ; '; else precmd=''; end

this_path = fileparts(which(mfilename));

% check if we use the cif2hkl executable, or need to compile the MeX (only once)
if ~isdeployed && isempty(compiled)
  compiled = '';
  
  % get list of modules
  modules = fullfile(this_path,'CFML*.f90');
  
  % check if iLoad config allows MeX
  config = iLoad('config');
  if isfield(config, 'MeX') && ...
   ((isfield(config.MeX, mfilename) && strcmp(config.MeX.(mfilename), 'yes')) ...
     || strcmp(config.MeX, 'yes'))
    % attempt to compile MeX
    fprintf(1, '%s: compiling mex...\n', mfilename);
    try
      cmd={'-c', '-O', '-output',fullfile(this_path,'cif2hkl.o'), ...
           modules, fullfile(this_path,'cif2hkl.F90')};
      disp([ 'mex ' sprintf('%s ', cmd{:}) ]);
      mex (cmd{:});
      cmd={'-O', '-output', fullfile(this_path,mfilename), ...
        fullfile(this_path,'cif2hkl_mex.c'), fullfile(this_path,'cif2hkl.o'), ...
        fullfile(this_path,'CFML*.o'), ...
        '-lgfortran'};
      disp([ 'mex ' sprintf('%s ', cmd{:}) ]);
      mex (cmd{:});
      compiled = 'mex';
    catch
      warning('%s: Can''t compile cif2hkl.F90 as MeX\n       in %s\n', ...
        mfilename, fullfile(this_path));
    end
    delete('*.mod');  % remove temprary files
    delete('*.o');
    if ~isempty(compiled)
      rehash
      % rethrow cif2hkl command with new MeX
      if exist(mfilename) == 3
        result = cif2hkl(varargin{:})
      elseif nargin > 0
        error('%s: MeX compiled successfully. Rethrow command to use it.');
      end
    end
    return
  else
    % look for executable, global and local
    [st, res] = system([ precmd mfilename ]);
    if st == 0
      compiled = mfilename; % exists globally in system PATH
    else
      % try local execution
      cmd = fullfile(this_path, mfilename);
      [st, res] = system([ precmd cmd ]);
      if st == 0
        compiled = cmd; % exists globally in system PATH
      else
        % attempt to compile as local binary
        if isempty(dir(fullfile(this_path,mfilename))) % no executable available
          fprintf(1, '%s: compiling binary...\n', mfilename);
          cmd = {'gfortran', '-O2', '-o', fullfile(this_path,mfilename), ...
             fullfile(this_path, 'CFML*.f90'), ...
             fullfile(this_path,'cif2hkl.F90'), '-lm'}; 
          disp([ sprintf('%s ', cmd{:}) ]);
          [status, result] = system(sprintf('%s ', cmd{:}));
          if status ~= 0 % not OK, compilation failed
            warning('%s: Can''t compile cif2hkl.F90 as binary\n       in %s\n', ...
              mfilename, fullfile(this_path));
          else
            compiled = fullfile(this_path,mfilename);
          end
        end
      end
    end
  end
end

if isempty(varargin) || strcmp(varargin{1}, 'compile'), return; end

% handle input arguments
% varargin = file_in, file_out, lambda, mode, verbose

if nargin < 1
  file_in = [];
else
  file_in = varargin{1};
end
if isempty(file_in) || ~ischar(file_in)
  error('At least one string input required. Syntax: cif2hkl(file_in, file_out, lambda, mode, verbose)')
end
if nargin < 2
  file_out = [];
else
  file_out = varargin{2};
end
if isempty(file_out)
  file_out = [ file_in '.laz' ];
end
if nargin < 3
  lambda = [];
else
  lambda = varargin{3};
end
if isempty(lambda) || ~isnumeric(lambda)
  lambda=0.5;
end
if nargin < 4
  mode = [];
else
  mode = varargin{4};
end
if isempty(mode) || ~ischar(mode)
  mode='p'; % powder or xtal
end
if nargin < 5
  verbose = [];
else
  verbose = varargin{5};
end
if isempty(verbose) || ~isnumeric(verbose)
  verbose=0;
end

% assemble command line
if ~isempty(compiled)
  cmd = compiled;
else
  cmd = mfilename;
end
if ~any(strcmp(file_in, {'-h','--help','--version'}))
  
  if verbose
    cmd = [ cmd ' --verbose' ];
  end
   % cmd = [ cmd ' --lambda ' num2str(lambda) ];
  if mode(1) == '-'
    cmd = [ cmd ' --no-outout-files' ];
  else
    if mode(1) == 'p'
      cmd = [ cmd ' --powder' ];
    elseif mode(1) == 'x'
      cmd = [ cmd ' --xtal' ];
    end
    cmd = [ cmd ' --out ' file_out ];
  end
  
end

cmd = [ cmd ' ' file_in ];
disp(cmd)

% launch the command
[status, result] = system(cmd);
if status ~= 0
  disp(result)
  error([ mfilename ' executable not available. Compile it with: "gfortran -O2 -o cif2hkl CFML_*.f90 cif2hkl.F90 -lm" in ' fullfile(fileparts(which(mfilename))) ]);
end


