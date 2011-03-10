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
%         dim: dimension to operate (int//array of)
% output: s: product of elements (iData/scalar)
% ex:     c=prod(a);
%
% Version: $Revision: 1.8 $
% See also iData, iData/plus, iData/prod, iData/cumprod, iData/mean

if ~isa(a, 'iData')
  iData_private_error(mfilename,['syntax is prod(iData, dim)']);
end

if nargin < 2, dim=1; end
% handle input iData arrays
if length(a(:)) > 1
  if dim ~= 0, s = a; else s=zeros(size(a)); end
  for index=1:length(a(:))
    s(index) = prod(a(index), dim);
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
% make axes single vectors for prod/trapz/... to work
for index=1:ndims(a)
  [x, xlab] = getaxis(a, index);
  setaxis(a, index, unique(x), xlab);
end

s = iData_private_cleannaninf(get(a,'Signal'));
e = iData_private_cleannaninf(get(a,'Error'));
m = iData_private_cleannaninf(get(a,'Monitor'));

[link, label] = getalias(a, 'Signal');
cmd= a.Command;
b  = copyobj(a);
rmaxis(b);

if all(dim > 0)
  % prod on all dimensions requested
  for index=1:length(dim(:))
    s = prod(s, dim(index)); 
    if numel(e) > 1, e = prod(e, dim(index)); end
    if numel(m) > 1, m = prod(m, dim(index)); end
  end
  % Store Signal
  s=squeeze(s); e=squeeze(e); m=squeeze(m);
  setalias(b,'Signal', s, [mfilename ' of ' label ]);
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

