function s = sum(a,dim)
% s = sum(a,dim) : computes the sum/projection of iData objects elements
%
%   @iData/sum function to compute the sum of the elements of the data set
%     sum(a,dim) accumulates along axis of rank dim. The axis is then removed.
%       If dim=0, sum is done on all axes and the total is returned as a scalar value. 
%       sum(a,1) accumulates on first dimension (columns)
%     sum(a,-dim) accumulates on all axes except the dimension specified, i.e.
%       the result is the projection of a along dimension dim.
%       All other axes are removed.
%
% input:  a: object or array (iData/array of)
%         dim: dimension to accumulate (int)
% output: s: sum of elements (iData/scalar)
% ex:     c=sum(a);
%
% See also iData, iData/plus, iData/prod, iData/cumsum, iData/mean

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
[link, label]          = getalias(a, 'Signal');
b=copyobj(a);
setaxis(b, [], getaxis(b)); % delete all axes
if all(dim > 0)
  for index=1:length(dim(:))
    s = sum(s, dim(index));
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
    s = sum(s, index);
  end
  return  % scalar
else  % dim < 0
  % accumulates on all axes except the rank specified
  for index=1:ndims(a)
    if index~=-dim, s = sum(s,index); end
  end
  setaxis(b, 1, getaxis(a, num2str(-dim)));
  setalias(b,'Signal', s, [ 'projection of ' label ]);     % Store Signal
end

b = iData_private_history(b, mfilename, b, dim);
s = b;

