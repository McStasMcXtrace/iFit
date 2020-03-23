function v = private_verbose(value)
% PRIVATE_VERBOSE set/get persistent/shared verbosity level
%   verbosity level
%     0 silent
%     1 normal
%     2 verbose
%     3 debug

  persistent verbose
  
  if isempty(verbose), verbose = 1; end
  
  if nargin
    verbose = value;
  end
  v = verbose;
