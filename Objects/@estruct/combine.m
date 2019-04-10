function c = combine(a,varargin)
% c = combine(a,b) : combines estruct objects
%
%   @estruct/combine (\) function to combine data sets
%     A fast notation for combine(a,b) is a\b
%     To combine a set of estruct objects use combine([ a b c ...])
%
%     The combine operator (aka merge) is defined as:
%       (S1+S2) over monitor (M1+M2)
%     where S1,M1 and S2,M2 and the Signal and Monitor of the two objects.
%
% input:  a: object or numerical array (estruct or numeric)
%         b: object or numerical array (estruct or numeric)
% output: c: object (estruct)
% ex:     c=combine(a,b); or combine([ a b ])
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


