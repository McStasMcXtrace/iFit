function s = trapz(a,dim)
% s = trapz(a,dim) : computes the integral of iData objects elements along given dimension
%
%   @iData/trapz function to compute the integral of the data set along a given dimension
%     trapz(a,dim) integrates along axis of rank dim. The axis is then removed.
%       default is to use dim=1
%       trapz is complementary to sum and camproj, but takes into account axis.
%
% input:  a: object or array (iData/array of)
%         dim: dimension to project (int)
% output: s: integral of elements (iData/scalar)
% ex:     c=trapz(a);
%
% Version: $Revision: 1.1 $
% See also iData, iData/cumsum, iData/camproj, iData/sum

if ~isa(a, 'iData')
  iData_private_error(mfilename,['syntax is trapz(iData, dim)']);
end

if nargin < 2, dim=1; end
% handle input iData arrays
if length(a(:)) > 1
  s = a;
  for index=1:length(a(:))
    s(index) = trapz(a(index), dim);
  end
  s = reshape(s, size(a));
  return
end

% we have a single object
s=get(a,'Signal');
[link, label]          = getalias(a, 'Signal');
[x, xlab]=getaxis(a,dim);
cmd=a.Command;
b=copyobj(a);
setaxis(b, [], getaxis(b)); % delete all axes

if all(dim > 0)
  for index=1:length(dim(:))
    if dim(index) ~= 1  % we put the dimension to integrate on as first
      perm=1:ndims(a);
      perm(dim(index))=1; perm(1)=dim(index);
      s = permute(s,perm);
    end
    s = trapz(x, s);
    if dim(index) ~= 1  % restore initial axes
      s = permute(s,perm);
    end
  end
  ax_index=1;
  for index=1:ndims(a)
    if all(index ~= dim)
      setaxis(b, ax_index, getaxis(a, num2str(index)));
      ax_index = ax_index+1;
    end
  end
  setalias(b,'Signal', s, [mfilename ' of ' label ' along ' xlab ]);     % Store Signal
elseif dim == 0
  for index=1:ndims(a)
    a = trapz(a, index);
  end
  return  % scalar
end
b.Command=cmd;
b = iData_private_history(b, mfilename, b, dim);
s = b;

