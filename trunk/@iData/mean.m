function a = mean(a)
% b = mean(s) : mean value of iData object
%
%   @iData/mean function to compute the mean value of objects
%
% input:  s: object or array (iData/array of)
% output: b: object or array (iData)
% ex:     b=mean(a);
%
% See also iData, iData/floor, iData/ceil, iData/round, iData/combine

if length(a) > 1
  a = combine(a);
end

b = mean(a.Signal);

