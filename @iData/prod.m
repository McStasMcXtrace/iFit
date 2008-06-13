function s = prod(a,dim)
% s = prod(a,dim) : computes the product of iData objects elements
%
%   @iData/prod function to compute the product of the elements of the data set
%     prod(a,dim) operates along axis of rank dim. The axis is then removed.
%       If dim=0, product is done on all axes and the total is returned as a scalar value. 
%       prod(a,1) operates on first dimension (columns)
%     prod(a,-dim) operates on all axes except the dimension specified, i.e.
%       the result is the product projection of a along dimension dim. 
%       All other axes are removed.
%
% input:  a: object or array (iData/array of)
%         dim: dimension to operate (int)
% output: s: product of elements (iData/scalar)
% ex:     c=prod(a);
%
% Version: $Revision: 1.5 $
% See also iData, iData/plus, iData/sum, iData/cumprod, iData/mean

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
[link, label]          = getalias(a, 'Signal');
b=copyobj(a);
setaxis(b, [], getaxis(b)); % delete all axes
if all(dim > 0)
  for index=1:length(dim(:))
    s = prod(s, dim(index));
  end
  ax_index=1;
  for index=1:ndims(a)
    if all(index ~= dim)
      setaxis(b, ax_index, getaxis(a, num2str(index)));
      ax_index = ax_index+1;
    end
  end
  setalias(b,'Signal', s, [mfilename ' of ' label ]);     % Store Signal
elseif dim == 0
  for index=1:ndims(a)
    s = prod(s, index);
  end
  return  % scalar
else  % dim < 0
  % operates on all axes except the rank specified
  for index=1:ndims(a)
    if index~=-dim, s = prod(s,index); end
  end
  setaxis(b, 1, getaxis(a, num2str(-dim)));
  setalias(b,'Signal', s, [ 'product projection of ' label ]);     % Store Signal
end

b = iData_private_history(b, mfilename, b, dim);
s = b;

