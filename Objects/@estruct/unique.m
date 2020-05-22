function s = unique(a,dim,mode)
% UNIQUE Set unique axes values with no repetitions.
%   S = UNIQUE(A) removes duplicates and sorts object first axis. The axis 
%   values are also sorted in ascending order.
%   Alternatively, you can use ISEQUAL to find unique data sets in a vector of
%   objects, and remove duplicates.
%
%   S = UNIQUE(A,DIM) sets unique values along axis of rank DIM. 
%   When DIM=0, operates on all axes.
%
% Example: a=estruct(1:21); a{1}=[9 9 9 9 9 9 8 8 8 8 7 7 7 6 6 6 5 5 4 2 1]; ...
%          b=unique(a); size(b,2) == 8
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plus, estruct/unique, estruct/sort, estruct/isequal

if nargin < 2, dim=1; end

% handle input estruct arrays
if numel(a) > 1
  s = zeros(estruct, numel(a), 1);
  for index=1:numel(a)
    s(index) = unique(a(index), dim);
  end
  s = reshape(s, size(a));
  return
end
cmd=a.Command;
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

% apply unique along axes
for index=tounique
  x = getaxis(a, index);
  [x, uniquei] = unique(x);
  if length(uniquei) ~= size(a, index)
    toeval='';
    S.type='()';
    for j=1:ndims(a), 
      if j ~= index, S.subs{j}=':';
      else           S.subs{j}=uniquei; end
    end
    % modify Signal, Error and Monitor
    if numel(se) == numel(sd), se =subsref(se, S); end
    if numel(sm) == numel(sd), sm =subsref(sm, S); end
    sd =subsref(sd, S);
    setaxis(s, dim, x);
    was_uniqueed=1;
  end
end
if was_uniqueed
  s = set(s, 'Signal', sd); label(s,'Signal',[ 'unique(' sl ')' ]);
  s = set(s, 'Error',  se);
  s = set(s, 'Monitor',sm);
  s.Command=cmd;
  s = history(s, mfilename, a, dim);
end

