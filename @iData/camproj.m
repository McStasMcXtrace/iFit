function s = camproj(a,dim)
% s = camproj(a,dim) : computes the projection of iData objects elements
%
%   @iData/camproj function to compute the projection/sum of the elements of the data set
%     camproj(a,dim) projects along axis of rank dim. All other axes are removed.
%       If dim=0, projection is done on all axes and the total is returned as a scalar value. 
%       camproj(a,1) projects on first dimension (rows).
%       camproj is the complementary to sum.
%
% input:  a: object or array (iData/array of)
%         dim: dimension to project (int)
% output: s: projection of elements (iData/scalar)
% ex:     c=camproj(a);
%
% Version: $Revision: 1.10 $
% See also iData, iData/plus, iData/prod, iData/cumsum, iData/mean, iData/sum, iData/trapz

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
  [x, xlab] = getaxis(a, index);
  setaxis(a, index, unique(x), xlab);
end

s = get(a,'Signal');
e = get(a,'Error');
m = get(a,'Monitor');

[link, label] = getalias(a, 'Signal');
cmd = a.Command;
b   = copyobj(a);
rmaxis(b, []); % removes all axes

if dim == 0
  for index=1:ndims(a)
    s = sum(s, index);
  end
  return  % scalar
else
  % accumulates on all axes except the rank specified
  for index=1:ndims(a)
    if index~=dim, 
      s = sum(s, index); 
      if numel(e) > 1, e = sum(e, index); e = sqrt(e.*e); end
      if numel(m) > 1, m = sum(m, index); end
    end
  end
  setalias(b,'Signal', s, [ 'projection of ' label ]);     % Store Signal
  b = set(b, 'Error', abs(e), 'Monitor', m);
  % set projection axis
  [x, xlab] = getaxis(a, num2str(dim)); % get axis definition and label
  setaxis(b, 1, x);
end

if dim == 1, 
	s=size(b);
	if s(1) == 1, b=transpose(b); end
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

