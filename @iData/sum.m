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
% Version: $Revision: 1.16 $
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
cmd= a.Command;
b  = copyobj(a);
rmaxis(b);

if all(dim > 0)
  % sum on all dimensions requested
  for index=1:length(dim(:))
    s = sum(s, dim(index)); 
    if numel(e) > 1, e = sum(e, dim(index)); e = sqrt(e.*e); end
    if numel(m) > 1, m = sum(m, dim(index)); end
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
  for index=1:ndims(a)
    s = feval(mfilename, a, index);
  end
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

