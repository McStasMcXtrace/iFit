function y = subsindex(s)
% d = subsindex(s) : subscript index for iData objects
%
%   @iData/subsindex: subscript index for iData objects
%   I = subsindex(S) is called for the syntax 'X(S)' to use S as an index.
%   The returned index is the convertion of S to integer.
%
% See also iData, iData/subsasgn, iData/subsref

% EF 23/09/07 iData implementation

if length(s) > 1
   s = s(1);
end

y = get(s,'Signal');
u = unique(y);
if all(u == 0 | u == 1)
  % signal is logical already
  y = find(y);
end
