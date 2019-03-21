function v=ndims(s)
% NDIMS get the dimensionality of estruct object (Signal).
%
% d = NDIMS(s) get the dimensionality of estruct object, i.e. that of its Signal.
% For an input array of objects, ndims returns the number of dimensions of the array.
%
% Example: ndims(estruct(1:10)) == 1 & ndims(estruct(peaks)) == 2
% Version: $Date$ (c) E.Farhi. License: EUPL.
% See also size

if numel(s) > 1
  v = arrayfun('ndims',s);
  return
end

n = size(s);
if     all(n == 0), v=0;
elseif all(n == 1), v=1;
else
  v = length(find(n > 1));
end


