function s = trapz(a,dim)
% s = trapz(a,dim) : computes the integral of iData objects elements along given dimension
%
%   @iData/trapz function to compute the integral of the data set along a given dimension
%     trapz(a,dim) integrates along axis of rank dim. The axis is then removed.
%       default is to use dim=1. If dim=0, integration is done on all axes and 
%       the total is returned as a scalar value. 
%       trapz is complementary to sum and camproj, but takes into account axis.
%
% input:  a: object or array (iData/array of)
%         dim: dimension to project (int)
% output: s: integral of elements (iData/scalar)
% ex:     c=trapz(a);
%
% Version: $Revision: 1.3 $
% See also iData, iData/cumsum, iData/camproj, iData/sum

if ~isa(a, 'iData')
  iData_private_error(mfilename,[ 'syntax is ' mfilename '(iData, dim)' ]);
end

if nargin < 2, dim=1; end
% handle input iData arrays
if length(a(:)) > 1
  s = a;
  for index=1:length(a(:))
    s(index) = feval(mfilename, a(index), dim);
  end
  s = reshape(s, size(a));
  return
end

% in all cases, resample the data set on a grid
a = interp(a,'grid');
% make axes single vectors for sum/trapz/... to work
for index=1:ndims(a)
  [x, lab] = getaxis(a, index);
  setaxis(a, index, unique(x),lab);
end

s = get(a,'Signal');
e = get(a,'Error');
m = get(a,'Monitor');

[link, label] = getalias(a, 'Signal');
cmd= a.Command;
b  = copyobj(a);
rmaxis(b, dim); % delete all axes to integrate, but keep any predefined alias

if all(dim > 0)
  for index=1:length(dim(:))
    [x, xlab]     = getaxis(a,dim(index));
    if dim(index) ~= 1  % we put the dimension to integrate on as first
      perm=1:ndims(a);
      perm(dim(index))=1; perm(1)=dim(index);
      s = permute(s, perm);
      e = permute(e, perm); 
      m = permute(m, perm);
    end
    % make the integration
    s = trapz(x, s);
    if numel(e) > 1, e = trapz(x, e); e = sqrt(e.*e); end
    if numel(m) > 1, m = trapz(x, m); end
    if dim(index) ~= 1  % restore initial axes
      s = permute(s,perm);
      e = permute(e,perm);
      m = permute(m,perm);
    end
  end
  % Store Signal
  setalias(b,'Signal', s, [mfilename ' of ' label ' along ' xlab ]);     
  b = set(b, 'Error', abs(e), 'Monitor', m);
elseif dim == 0
  for index=1:ndims(a)
    s = feval(mfilename, a, index);
  end
  return  % scalar
end
b.Command=cmd;
b = iData_private_history(b, mfilename, b, dim);
s = b;

