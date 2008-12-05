function s = sort(a,dim,mode)
% s = sort(a,dim) : Sort iData objects axes in ascending or descending order
%
%   @iData/sort function to sort the data set on its axes
%     sort(a,dim) sorts along axis of rank dim. 
%       If dim=0, sorting is done on all axes.
%     sort(a,dim,mode) where mode='ascend' or 'descend' select sorting order
%
% input:  a: object or array (iData)
%         dim: dimension to sort (int)
%         mode: sorting order 'ascend' or 'descend'
% output: s: sorted data (iData)
% ex:     c=sort(a);
%
% Version: $Revision: 1.4 $
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
se = get(s,'Error');
sm = get(s,'Monitor');

if dim > 0
  tosort=dim;
else
  tosort=1:ndims(a)
end
was_sorted=0;

for index=tosort
  x = getaxis(a, index);
  [x, sorti] = sort(x, index, mode);
  if ~isequal(sorti, 1:size(a, index))
    toeval='';
    for j=1:ndims(a), 
      if j ~= index, str_idx{j}=':';
      else str_idx{j}='sorti'; end
      if j>1, toeval=[ toeval ',' str_idx{j} ];
      else toeval=[ str_idx{j} ]; end
    end
    try
      sd =eval([ 'sd(' toeval ')' ]);
    catch
    end
    try
      se =eval([ 'se(' toeval ')' ]);
    catch
      sd=[];
    end
    try
      sm =eval([ 'sm(' toeval ')' ]);
    catch
      sm=[];
    end
    setaxis(s, index, x);
    was_sorted=1;
  end
end
if was_sorted
  s = setalias(s, 'Signal', sd, [ 'sort(' sl ')' ]);
  s = setalias(s, 'Error',  se);
  s = setalias(s, 'Monitor',sm);
  s = iData_private_history(s, mfilename, a, dim, mode);
end

