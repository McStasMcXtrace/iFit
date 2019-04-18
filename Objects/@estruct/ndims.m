function v=ndims(s)
% NDIMS get the dimensionality of estruct object (Signal).
%   D = NDIMS(s) get the dimensionality of object, i.e. that of its Signal.
%   For an input array of objects, ndims returns the number of dimensions
%   of the array.
%
%   To get the dimensionality of all objects in an array, use ARRAYFUN('ndims',S)
%
% Example: ndims(estruct(1:10)) == 1 & ndims(estruct(peaks)) == 2
% Version: $Date$ $Version$ $Author$
% See also size

if numel(s) > 1  % this is an array of estruct
  v = builtin('size', s);
  return
end

n = size(s);
if     all(n == 0), v=0;
elseif all(n == 1), v=1;
else
  v = length(find(n > 1));
end