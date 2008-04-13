function s = prod(a,dim)
% s = prod(a,dim) : computes the prod of iData objects elements
%
%   @iData/prod function to compute the product of the elements of the data set
%     prod(a,dim) operates along axis of rank dim. If dim=0, product is done
%       on all axes.
%
% input:  a: object or array (iData)
%         dim: dimension to operate (int)
% output: s: prod of elements (double)
% ex:     c=prod(a);
%
% See also iData, iData/plus, iData/sum, iData/cumprod

if ~isa(a, 'iData')
  iData_private_error(mfilename,['syntax is prod(iData, dim)']);
end

if nargin < 2, dim=1; end
% handle input iData arrays
if length(a(:)) > 1
  s = {};
  for index=1:length(a(:))
    s{index} = prod(a(index), dim);
  end
  s = reshape(s, size(a));
  return
end

s=get(a,'Signal');
if dim > 0
  s = prod(s, dim);
else
  for index=1:ndims(a)
    s = prod(s, index);
  end
end

