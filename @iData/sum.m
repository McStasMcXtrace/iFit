function s = sum(a,dim)
% s = sum(a,dim) : computes the sum of iData objects elements
%
%   @iData/sum function to compute the sum of the elements of the data set
%     sum(a,dim) accumulates along axis of rank dim. The axis is then removed.
%       If dim=0, sum is done on all axes and the total is returned as a scalar value. 
%       sum(a,1) accumulates on first dimension (columns). 
%       camproj accumulates on all other axes.
%
% input:  a: object or array (iData/array of)
%         dim: dimension to accumulate (int)
% output: s: sum of elements (iData/scalar)
% ex:     c=sum(a);
%
% Version: $Revision: 1.11 $
% See also iData, iData/plus, iData/prod, iData/cumsum, iData/mean, iData/camproj, iData/trapz

if ~isa(a, 'iData')
  iData_private_error(mfilename,['syntax is sum(iData, dim)']);
end

if nargin < 2, dim=1; end
% handle input iData arrays
if length(a(:)) > 1
  s = a;
  for index=1:length(a(:))
    s(index) = sum(a(index), dim);
  end
  s = reshape(s, size(a));
  return
end

s=get(a,'Signal');
[link, label]          = getalias(a, 'Signal');
cmd=a.Command;
b=copyobj(a);
setaxis(b, [], getaxis(b)); % delete all axes
if all(dim > 0)
  for index=1:length(dim(:))
    s = trapz(s, dim(index));
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
    s = trapz(s, index);
  end
  return  % scalar
end
b.Command=cmd;
b = iData_private_history(b, mfilename, b, dim);
s = b;

