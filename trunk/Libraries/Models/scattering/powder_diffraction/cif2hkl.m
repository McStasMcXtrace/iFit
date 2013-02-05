function result = cif2hkl(file_in, file_out, lambda, mode, verbose)
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

persistent compiled

this_path = fileparts(which(mfilename));

% handle input arguments
if nargin < 1
  file_in = [];
end
if isempty(file_in) || ~ischar(file_in)
  error('At least one string input required. Syntax: cif2hkl(file_in, file_out, lambda, mode, verbose)')
end
if nargin < 2
  file_out = [];
end
if isempty(file_out)
  file_out = [ file_in '.laz' ];
end
if nargin < 3
  lambda = [];
end
if isempty(lambda) || ~isnumeric(lambda)
  lambda=0.5;
end
if nargin < 4
  mode = [];
end
if isempty(mode) || ~ischar(mode)
  mode='p'; % powder or xtal
end
if nargin < 5
  verbose = [];
end
if isempty(verbose) || ~isnumeric(verbose)
  verbose=0;
end

% check if we use the cif2hkl executable, or need to compile the MeX
if ~isdeployed && isempty(compiled)
  p = pwd;
  cd (fullfile(fileparts(which(mfilename))))
  fprintf(1, '%s: compiling cif2hkl...\n', mfilename);
  try
    mex -c -O cif2hkl.F90
    mex -O cif2hkl_mex.c cif2hkl.o -o cif2hkl -lgfortran
    compiled = 1;
  catch
    compiled = 0;
    error('%s: Can''t compile cif2hkl.F90\n       in %s\n', ...
      mfilename, pwd);
  end
  cd (p)
  if compiled
    rehash
    error('%s: Compilation went FINE. Please rethrow your command so that the new MeX can be used...', mfilename);
  end
end

% assemble command line
cmd = fullfile(this_path, 'cif2hkl');
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
  error([ 'cif2hkl executable not available. Compile it with: "gfortran -O2 -o cif2hkl cif2hkl.f90 -lm" in ' fullfile(fileparts(which(mfilename))) ]);
end


