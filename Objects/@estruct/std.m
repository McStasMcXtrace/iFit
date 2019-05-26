function [w,c] = std(a, dim)
% STD Standard deviation of object.
%  [W, C] = STD(A) computes the standard deviation of object A, that is
%  their gaussian half width W (second moment). Optionally, the distribution  
%  center C (first moment) can be returned as well.
%  The minimum value object value (background) is subtracted before estimating
%  moments.
%
%  [W,C] = STD(A, DIM) computes standard deviation along axis of rank 'dim'.
%  When omitted, dim is set to 1. DIM can be given as an array to iteratively
%  compute the width and center of projected distributions along all given axes.
%
%  [W,C] = STD(A, -DIM) does not subtract the minimum object value (background) before 
%  computing the width and center. This may result in imaginary values.
%
%  [W,C] = STD(A, 0) computes std(Signal) and mean(Signal).
%
% Example: a=estruct(peaks); c=std(a);
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/median, estruct/mean
w = []; c = [];
if nargin < 2, dim=1; end
if numel(a) > 1
  for index=1:numel(a)
    [wi, ci] = std(a(index), dim);
    w = [ w wi ];
    c = [ c ci ];
  end
  return
end

if length(dim) > 1
  for index=1:length(dim)
    [wi, ci] = std(a, dim(index));
    w = [ w wi ];
    c = [ c ci ];
  end
  return
end

if abs(dim) > ndims(a)
  dim = 1;
end

if dim == 0
  w = double(a);
  c = mean(w(:));
  w = std(w(:));
  return
end

% use cache when available for faster execution
if isfield(a.Private,'cache') && dim > 0 ...
  if isfield(a.Private.cache,'std_w') && ~isempty(a.Private.cache.std_w)
    w = a.Private.cache.std_w;
    if numel(w) >= dim && w(dim), w = w(dim); else w = []; end
  end
  if  isfield(a.Private.cache,'std_c') && ~isempty(a.Private.cache.std_c)
    c = a.Private.cache.std_c;
    if numel(c) >= dim && c(dim), c = c(dim); else c = []; end
  end
end

if ~isempty(w) && (nargout == 1 || ~isempty(c))
  return
end

% we first compute projection of estruct on the selected dimension
if ~isvector(a) || length(a.Axes) < ndims(a)
  b = meshgrid(a); % make it a clean grid style data set
else b = a;
end

s = getaxis(b,'Signal'); 
x = getaxis(b, abs(dim));

s=s(:); x=x(:);
if ~isfloat(s), s=double(s); end
if ~isfloat(x), s=double(x); end

% then we compute sum(axis{dim}.*Signal)/sum(Signal)
s = private_cleannaninf(s);
if (dim > 0)
  s = s - min(s);
end

sum_s = sum(s);
% first moment (mean)
c = sum(s.*x)/sum_s; % mean value

% second moment: sqrt(sum(x^2*s)/sum(s)-fmon_x*fmon_x);
w = sqrt(sum(x.*x.*s)/sum_s - c*c);

if dim > 0
  % store in cache
  a.Private.cache.std_w(dim) = w;
  a.Private.cache.std_c(dim) = c;
elseif ~isreal(w) && dim < 0
  inputname1 = inputname(1);
  if isempty(inputname1), inputname1 = class(a); end

  warning([ mfilename ': The computed standard deviation is imaginary. ' ...
    'You should use std(' inputname1 ', ' num2str(-dim) ')' ])
end

