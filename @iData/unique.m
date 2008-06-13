function s = unique(a,dim,mode)
% s = unique(a,dim) : set unique iData objects axes with no repetitions
%
%   @iData/unique function to set unique the data set on its axes
%     unique(a,dim) set unique along axis of rank dim. 
%       If dim=0, operates on all axes.
%
% input:  a: object or array (iData)
%         dim: dimension to unique (int)
% output: s: data set with unique axes (iData)
% ex:     c=unique(a);
%
% Version: $Revision: 1.2 $
% See also iData, iData/plus, iData/unique, iData/sort
if ~isa(a, 'iData')
  iData_private_error(mfilename,['syntax is unique(iData, dim)']);
end

if nargin < 2, dim=1; end

% handle input iData arrays
if length(a(:)) > 1
  s = a(:);
  for index=1:length(a(:))
    s(index) = unique(a(index), dim);
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
  tounique=dim;
else
  tounique=1:ndims(a)
end
was_uniqueed=0;

for index=tounique
  x = getaxis(a, index);
  [x, uniquei] = unique(x);
  if length(uniquei) ~= size(a, index)
    toeval='';
    for j=1:ndims(a), 
      if j ~= index, str_idx{j}=':';
      else str_idx{j}='uniquei'; end
      if j>1, toeval=[ toeval ',' str_idx{j} ];
      else toeval=[ str_idx{j} ]; end
    end
    sd =eval([ 'sd(' toeval ')' ]);
    se =eval([ 'se(' toeval ')' ]);
    sm =eval([ 'sm(' toeval ')' ]);
    setaxis(s, dim, x);
    was_uniqueed=1;
  end
end
if was_uniqueed
  s = setalias(s, 'Signal', sd, [ 'unique(' sl ')' ]);
  s = setalias(s, 'Error',  se);
  s = setalias(s, 'Monitor',sm);
  s = iData_private_history(s, mfilename, a, dim);
end

