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
% Version: $Revision: 1.12 $
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

% in all cases, resample the data set on a grid
a = interp(a,'grid');
% make axes single vectors for sum/trapz/... to work
for index=1:ndims(a)
  x = getaxis(a, index);
  setaxis(a, index, unique(x));
end

s = get(a,'Signal');
e = get(a,'Error');
m = get(a,'Monitor');

[link, label] = getalias(a, 'Signal');
cmd= a.Command;
b  = copyobj(a);
setaxis(b, [], getaxis(b)); % delete all axes, but keep any predefined alias
if all(dim > 0)
  % sum on all dimensions requested
  for index=1:length(dim(:))
    s = sum(s, dim(index)); 
    if numel(e) > 1, e = sum(e, dim(index)); e = sqrt(e.*e); end
    if numel(m) > 1, m = sum(m, dim(index)); end
  end
  % reconstruct all required axes, except the one removed
  ax_index=1;
  for index=1:ndims(a)
    if all(index ~= dim)  % copy all axes except those which are summed
      x = getaxis(a, index);
      setaxis(b, ax_index, x);
      ax_index = ax_index+1;
    end
  end
  % Store Signal
  setalias(b,'Signal', s, [mfilename ' of ' label ]);
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

