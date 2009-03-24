function c = combine(a,varargin)
% c = combine(a,b) : combines iData objects
%
%   @iData/combine (\) function to combine data sets
%     A fast notation for combine(a,b) is a\b
%     To combine a set of iData objects use combine([ a b c ...])
%
% input:  a: object or numerical array (iData or numeric)
%         b: object or numerical array (iData or numeric)
% output: c: object (iData)
% ex:     c=combine(a,b); or combine([ a b ])
%
% Version: $Revision: 1.6 $
% See also iData, iData/minus, iData/plus, iData/times, iData/rdivide
if length(varargin) > 1  % syntax: combine(a,b,...)
  s=a(:);
  for index=1:length(varargin)
    s = [ s ; varargin{index} ];
  end
  c = combine(s);
  return
end

% now we should only handle a single argument
c = a(1);
if length(a) <= 1, return; end
for index=2:length(a)
  c = iData_private_binary(c, a(index), 'combine');
end


