function c = colon(a,d,b)
% :  Colon (linear morphing).
%   C = J:K = COLON(J,K) is the same as [J, J+1, ..., K], i.e. creates a vector
%   of objects that goes from J to K, using LINSPACE.
%
% Example: b=iData(peaks); a=-b; c=a:.5:b; numel(c) == 4
% Version: $Date$ $Version$ $Author$
% See also iData, iData/floor, iData/ceil, iData/round, iData/combine

if nargin == 1
  b = []; d=0;
elseif nargin == 2
  b = d;
  d = 0;
end
if isempty(b)
  b=a;
end
if numel(a) > 1, a=a(1); end
if numel(b) > 1, b=b(end); end
if isempty(d), d=0; end

if isa(a, 'iData'), as = getaxis(a(1),'Signal'); as = as(:); else as = a; end
if isa(b, 'iData'), bs = getaxis(b(1),'Signal'); bs = bs(:); else bs = b; end

as = floor(mean(as));
bs = ceil( mean(bs));
if d    , n = abs((bs-as)/d);
else      n = abs( bs-as); end

if n <= 1, c = a; return; end
if n < 0,  c = a; return; end

c = linspace(a,b,n);

