function v = isvector(s)
% ISVECTOR True if array is a vector.
%   ISVECTOR(V) returns logical true (1) if V is a 1 x n or n x 1 vector,
%   where n >= 0, and logical false (0) otherwise.
%   Additionally, if the vector data set has D axes defined, it returns D, for
%   instance with event data sets.
%
% Example: s=estruct(1:10); isvector(s) == 1
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sign, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar,
%          estruct/isvector, estruct/issparse


v = unary(s, 'isvector');
if iscell(v), v = cell2mat(v); end

% special case for event data sets: signal is a vector, all axes have same dimension
for index=1:numel(s) % OK with object arrays
  if v(index) == 1 && numel(getaxis(s(index))) > 1

    % check that all axes have same length, else keep v=1
    l = prod(size(s(index)));
    flag = true;
    for dim=1:numel(getaxis(s(index)))
      if length(getaxis(s(index),dim)) ~= l, flag=false; break; end
    end % dim
    if flag, v(index) = 1+getaxis(s(index)); end % [Signal, Axes]

  end

end % index
