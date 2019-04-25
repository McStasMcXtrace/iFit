function c = combine(a,varargin)
% COMBINE Combines (merge) objects Signal and axes.
%   C = COMBINE(A,B) combines all objects in [A B] into a single object.
%   The merge is carried out on Signal and Axes.
%   The combine operator (aka merge) of two objects is defined as:
%     (S1+S2) over monitor (M1+M2)
%   where S1,M1 and S2,M2 and the Signal and Monitor of the two objects.
%
%   Alternatively, the addition (PLUS) is defined as the normalised sum:
%       (M1+M2)*(S1/M1+S2/M2) over monitor(M1+M2)
%
%   C = COMBINE([ A B C ...]) combines iteratively all objects.
%   C = COMBINE(A, B, C, ...) combines iteratively all objects.
%
%   A\B and A.\B are equivalent to COMBINE(A,B).
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/minus, estruct/plus, estruct/times, estruct/rdivide

if length(varargin) >= 1  % syntax: combine(a,b,...)
  s=a(:);
  for index=1:length(varargin)
    s = [ s ; varargin{index} ];
  end
  clear varargin
  c = combine(s);
  return
end

% now we should only handle a single argument
if all(isvector(a) > 1)
  % combine event data sets by simple catenation
  c = cat(1, a);
else
  c = a(1);
  if length(a) <= 1, return; end
  c = binary(a, [], 'combine');
end