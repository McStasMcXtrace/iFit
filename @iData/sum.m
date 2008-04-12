function s = sum(a,dim)
% s = sum(a,dim) : computes the sum of iData objects elements
%
%   @iData/sum function to compute the sum of the elements of the data set
%     sum(a,dim) accumulates along axis of rank dim.
%
% input:  a: object or array (iData or numeric)
%         dim: dimension to accumulate (int)
% output: s: sum of elements (double)
% ex:     c=sum(a);
%
% See also iData, iData/plus, iData/prod, iData/cumsum

% handle input iData arrays
if isa(a, 'iData') & length(a(:)) > 1
  b = zeros(length(a(:)));
  for index=1:length(a(:))
    if nargin == 1
      b(index) = sum(a);
    else
      b(index) = sum(a, dim);
    end
  end
  b = reshape(b, size(a));
  return
end

if nargin == 1
  s = sum(get(a,'Signal'));
else
  if isa(a, 'iData')
    s = sum(get(a,'Signal'), dim);
  else
    iData_private_error(mfilename,['syntax is sum(iData, dim)']);
  end
end

