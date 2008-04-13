function s = sort(a,dim,mode)
% s = sort(a,dim) : Sort iData objects elements in ascending or descending order
%
%   @iData/sort function to compute the cumulative sum of the elements of the data set
%     sort(a,dim) sorts along axis of rank dim. If dim=0, sorting is done
%       on all axes.
%     sort(a,dim,mode) where mode='ascend' or 'descend' select sorting order
%
% input:  a: object or array (iData)
%         dim: dimension to sort (int)
%         mode: sorting order 'ascend' or 'descend'
% output: s: sorted data (iData)
% ex:     c=sort(a);
%
% See also iData, iData/plus, iData/sort, iData/unique
if ~isa(a, 'iData')
  iData_private_error(mfilename,['syntax is sort(iData, dim, mode)']);
end

if nargin < 2, dim=1; end
if nargin < 3, mode='ascend'; end

% handle input iData arrays
if length(a(:)) > 1
  s = a(:);
  for index=1:length(a(:))
    s(index) = sort(a(index), dim, mode);
  end
  s = reshape(s, size(a));
  return
end

s = copyobj(a);

[sn, sl] = getaxis(a, '0');   % label
sd = get(s,'Signal');         % data

if dim > 0
  [sd, sorti] = sort(sd, dim, mode);
  if ~isempty(getaxis(s,num2str(dim)))
    x = getaxis(s, dim);
    setaxis(s, dim, x(sorti));
  end
else
  for index=1:ndims(a)
    [sd, sorti] = sort(sd, index, mode);
    if ~isempty(getaxis(s,num2str(index)))
      x = getaxis(s, index);
      setaxis(s, index, x(sorti));
    end
  end
end

s = setalias(s, 'Signal', sd, [ 'sort(' sl ')' ]);
s = iData_private_history(s, mfilename, a, dim, mode); 

