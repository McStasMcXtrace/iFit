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
%         dim: dimension to integrate (int//array of)
% output: s: integral of elements (iData/scalar)
% ex:     c=trapz(a);
%
% Version: $Revision: 1.8 $
% See also iData, iData/cumsum, iData/camproj, iData/sum

if ~isa(a, 'iData')
  iData_private_error(mfilename,[ 'syntax is ' mfilename '(iData, dim)' ]);
end

if nargin < 2, dim=1; end
% handle input iData arrays
if length(a(:)) > 1
  if dim ~= 0, s = a; else s=zeros(size(a)); end
  for index=1:length(a(:))
    s(index) = feval(mfilename, a(index), dim);
  end
  s = reshape(s, size(a));
  return
end

% removes warnings
try
  warn.set = warning('off','iData:setaxis');
  warn.get = warning('off','iData:getaxis');
catch
  warn = warning('off');
end

% in all cases, resample the data set on a grid
a = interp(a,'grid');
% make axes single vectors for sum/trapz/... to work
for index=1:ndims(a)
  [x, lab] = getaxis(a, index);
  setaxis(a, index, unique(x),lab);
end

s = iData_private_cleannaninf(get(a,'Signal'));
e = iData_private_cleannaninf(get(a,'Error'));
m = iData_private_cleannaninf(get(a,'Monitor'));

[link, label] = getalias(a, 'Signal');
cmd= a.Command;
b  = copyobj(a);
rmaxis(b);

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
  % put back initial axes, except those integrated
  ax_index=1;
  for index=1:ndims(a)
    if all(dim ~= index)
      [x, xlab] = getaxis(a, num2str(index)); % get axis definition and label
      setaxis(b, ax_index, x);
      ax_index = ax_index+1;
    end
  end
elseif dim == 0
  s = feval(mfilename, a, 1:ndims(a));
  s = double(s);
  return  % scalar
end
b.Command=cmd;
b = iData_private_history(b, mfilename, b, dim);
s = b;

% reset warnings
try
  warning(warn.set);
  warning(warn.get);
catch
  warning(warn);
end

