function a = subsasgn_recursive(a, S, val, a0)
% subsasgn_recursive assign value recursively
%   SUBSASGN_RECURSIVE(a,S,val,a0) assigns recursively a.(S) = val in root object a0
  if nargin < 4, a0=a; end
  
  if numel(S) == 1 
    % set the value to the last level
    a = subsasgn_single(a, S, val, a0);
  else
    % will use our own subsref/subsasgn for iData, 
    % but default for other objects
    try
      b = subsref_single(a, S(1),a0);               % b         = a.(S(1)
      b = subsasgn_recursive(b, S(2:end), val, a0); % b.(S(2:)) = val
    catch
      % a0.( S(1) ) does not exist yet. We create it as empty
      if S(1).type == '.'
        f = S(1).subs;
        if iscellstr(f), f=char(f); end
        if isa(a, 'iData')
          addprop(a, f);
        else
          a.(f) = [];
        end
        b = [];
        b = subsasgn(b, S(2:end), val);
      end % silently ignore wrong assignments...
    end
    
    a = subsasgn_single(a, S(1), b, a0);
  end

