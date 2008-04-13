function b = mean(a, dim)
% b = mean(s) : mean value of iData object
%
%   @iData/mean function to compute the mean value of objects
%
% input:  s: object or array (iData/array of)
%         dim: dimension to average (int)
% output: b: object or array (iData)
% ex:     b=mean(a);
%
% See also iData, iData/floor, iData/ceil, iData/round, iData/combine

if nargin < 2, dim=1; end
if length(a) > 1
  a = combine(a);
  return
end

c = double(a);
b = mean(c(:), dim);

