function s = prod(a,dim)
% s = prod(a,dim) : computes the prod of iData objects elements
%
%   @iData/prod function to compute the product of the elements of the data set
%     prod(a,dim) operates along axis of rank dim.
%
% input:  a: object or array (iData or numeric)
%         dim: dimension to operate (int)
% output: s: prod of elements (double)
% ex:     c=prod(a);
%
% See also iData, iData/plus, iData/sum, iData/cumprod

% handle input iData arrays
if isa(a, 'iData') & length(a(:)) > 1
  b = zeros(length(a(:)));
  for index=1:length(a(:))
    if nargin == 1
      b(index) = prod(a);
    else
      b(index) = prod(a, dim);
    end
  end
  b = reshape(b, size(a));
  return
end

if nargin == 1
  s = prod(get(a,'Signal'));
else
  if isa(a, 'iData')
    s = prod(get(a,'Signal'), dim);
  else
    iData_private_error(mfilename,['syntax is prod(iData, dim)']);
  end
end

