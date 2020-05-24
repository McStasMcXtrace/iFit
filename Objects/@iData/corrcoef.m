function [r,p,rlo,rup] = corrcoef(a, b)
% CORRCOEF Correlation coefficients.
%   R = CORRCOEF(A,B) computes the correlation coefficient matrix between 
%   an object and other data. C ranges between -1 (opposite) and 1 (similar).
%   To get a single correlation measure, use the off-diagonal values e.g. R(1,2).
%   http://en.wikipedia.org/wiki/Correlation_coefficient
%
%   R = CORRCOEF(A) computes the correlation coefficient between an object
%   Signal and its Axes.
%
%   [R,P]=CORRCOEF(...) also returns P, a matrix of p-values for testing
%   the hypothesis of no correlation. Each p-value is the probability of getting
%   a correlation as large as the observed value by random chance, when the true
%   correlation is zero. If P(i,j) is small, say less than 0.05, then the 
%   correlation R(i,j) is significant.

%   [R,P,RLO,RUP]=CORRCOEF(...) also returns matrices RLO and RUP, of the same 
%   size as R, containing lower and upper bounds for a 95% confidence interval 
%   for each coefficient.
%
% Example: a=iData(peaks); r=corrcoef(a, -a); r(1,2) == -1
% Version: $Date$ $Version$ $Author$
% See also corrcoef, iData, iData/mean, iData/fits

% either 'a' or 'b' is an iData
if nargin > 1 && isa(b, 'iData')
  tmp = a; a = b; b= tmp;
end

% handle input iData arrays
if numel(a) > 1
  c = cell(numel(a),1);
  for index=1:numel(a)
    c{index} = feval(mfilename, a(index), b);
  end
  c = reshape(s, size(a));
  return
end

% handle one input case: use axis and signal
if nargin == 1
  if ndims(a) > 1
    a = meshgrid(a);
  end
  b = getaxis(a, 0);
  a = getaxis(a, 0);
end

% handle input iFunc arrays
if isa(b, 'iFunc')
  b = feval(b, NaN, a);
end

% find intersection
if isa(a, 'iData') && isa(b, 'iData')
  [a,b] = intersect(a,b); % perform operation on intersection
end
% get the Signal of the two objects
if isa(b, 'iData')
  b = getaxis(b, 'Signal');
end
if isa(a, 'iData')
  a = getaxis(a, 'Signal');
end

if ~isnumeric(a) || ~isnumeric(b) || numel(a) ~= numel(b)
  c = [];
  return
end

index = find(isfinite(a) & isfinite(b));
[r,p,rlo,rup] = corrcoef(a(index), b(index));

