function a = subsasgn_recursive(a, S, val, a0)
  % subsasgn_recursive assign value recursively
  if nargin < 4, a0=a; end
  
  if numel(S) == 1 
    % set the value to the last level
    a = subsasgn_single(a, S, val, a0);
  else
    % will use our own subsref/subsasgn for estruct, 
    % but default for other objects
    b = subsref_single(a, S(1),a0); 
    b = subsasgn_recursive(b, S(2:end), val, a0);
    a = subsasgn_single(a, S(1), b, a0);
  end

