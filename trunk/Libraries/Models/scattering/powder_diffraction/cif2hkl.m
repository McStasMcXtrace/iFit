function result = cif2hkl(file_in, file_out, lambda, mode, verbose)
% result = cif2hkl(file_in, file_out, lambda, mode, verbose)
%
% Input:
%   file_in:   name of CIF/CFL/PCR/ShellX file (char)
%   file_out:  name of output file, or empty for '<file_in>.laz' (char)
%   lambda:    minimum wavelength used for generation of HKL peaks. Default is 0.5  (scalar, Angs)
%   mode:      file conversion mode 'p'=powder, 'x'=single crystal, '-'=no/unactivated
%   verbose:   verbosity level, 0 or 1 (scalar)
% Output:
%   result:    command result (char)
%   file_out when mode is not '-'
%   

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
  error('cif2hkl executable not available. Compile it with: gfortran -O2 -o cif2hkl cif2hkl.f90 -lm');
end


