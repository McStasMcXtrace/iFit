function v = isvector(s)
% b = isvector(s) : True for vector iData object
%
%   @iData/isvector function to return true if data set is a vector
%   i.e. that its size is 1xN or Nx1.
%   Additionally, if the vector data set has X axes defined, it returns X.
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=isvector(a);
%
% See also iData, iData/sign, iData/isreal, iData/isfinite, iData/isnan,
%          iData/isinf, iData/isfloat, iData/isinterger,
%          iData/isnumeric, iData/islogical, iData/isscalar, 
%          iData/isvector, iData/issparse

if length(s(:)) > 1
  v = zeros(size(s)); 
  for index=1:length(s(:))
    v(index) =ndims(s(index));
  end
  return
end

n = size(s);
v = 0;
if all(n == 1), v=1;
else
  index=find(n > 1);
  v = length(index);
  if v == 1 & length(getaxis(s)) > 1
    v = length(getaxis(s)); % this is for [x,y,z,... vector data (plot3 style)]
  elseif v ~= 1, v=0; end
end

