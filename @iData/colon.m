function c = colon(a,d,b)
% b = mean(s) : mean value of iData object
%
%   @iData/mean function to compute the mean value of objects
%
% input:  s: object or array (iData/array of)
% output: b: object or array (iData)
% ex:     b=mean(a);
%
% Version: $Revision: 1.2 $
% See also iData, iData/floor, iData/ceil, iData/round, iData/combine

if nargin == 2
  b = d;
  d = 0;
end

as = round(mean(a));
bs = round(mean(b));
if d > 0, n = abs(bs-as)/d;
else      n = abs(bs-as); end

if n == 1, c = a;  return; end
if n <= 0, c = a; return; end

c = linspace(a,b,n);

