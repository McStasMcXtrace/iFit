function s = sum(a,dim)
% s = sum(a,dim) : computes the sum of iData objects elements
%
%   @iData/sum function to compute the sum of the elements of the data set
%     sum(a,dim) accumulates along axis of rank dim. If dim=0, sum is done
%       on all axes.
%
% input:  a: object or array (iData or numeric)
%         dim: dimension to accumulate (int)
% output: s: sum of elements (double)
% ex:     c=sum(a);
%
% See also iData, iData/plus, iData/prod, iData/cumsum

if ~isa(a, 'iData')
  iData_private_error(mfilename,['syntax is sum(iData, dim)']);
end

if nargin < 2, dim=1; end
% handle input iData arrays
if length(a(:)) > 1
  s = {};
  for index=1:length(a(:))
    s{index} = sum(a(index), dim);
  end
  s = reshape(s, size(a));
  return
end

s=get(a,'Signal');
if dim > 0
  s = sum(s, dim);
else
  for index=1:ndims(a)
    s = sum(s, index);
  end
end

