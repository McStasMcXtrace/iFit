function c = colon(a,d,b)
% b = colon(s) : vector of arrays
%
%   @iData/mean function to compute the mean value of objects
%
% input:  s: object or array (iData/array of)
% output: b: object or array (iData)
% ex:     b=mean(a);
%
% Version: $Revision: 1.3 $
% See also iData, iData/floor, iData/ceil, iData/round, iData/combine

if nargin == 2
  b = d;
  d = 0;
end

if isa(a, 'iData'), as = getaxis(a(1),0); as = as(:); else as = a; end
if isa(b, 'iData'), bs = getaxis(b(1),0); bs = bs(:); else bs = b; end

as = round(mean(as));
bs = round(mean(bs));
if d > 0, n = abs(bs-as)/d;
else      n = abs(bs-as); end

if n == 1, c = a;  return; end
if n <= 0, c = a; return; end

c = linspace(a,b,n);

