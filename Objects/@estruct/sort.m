function b = sort(a,dim,mode)
% SORT  Sort axes in ascending or descending order.
%   B = SORT(A) sorts first rank axes in ascending values. Object B is usually 
%   the same as A, but with ascending axes values.
%
%   B = SORT(A,DIM) sorts axis of rank DIM. 
%   When DIM=0, sorting is done sequentially on all axes.
%
%   B = SORT(A,DIM,MODE) where MODE='ascend' or 'descend' selects sorting order.
%
% Example: a=estruct(1:10); a{1}=-a{1}; b=sort(a); all(b{1}+a{1} == -11)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plus, estruct/sort, estruct/unique

if nargin < 2, dim=1; end
if nargin < 3, mode='ascend'; end

% handle input estruct arrays
if numel(a) > 1
  b = zeros(estruct, numel(a), 1);
  for index=1:numel(a)
    b(index) = sort(a(index), dim, mode);
  end
  b = reshape(b, size(a));
  return
end
cmd=a.Command;
b = copyobj(a);

sd = subsref(b,struct('type','.','subs','Signal'));
[dummy, sl] = getaxis(b, '0');  % signal definition/label
se = subsref(b,struct('type','.','subs','Error'));
sm = subsref(b,struct('type','.','subs','Monitor'));
if numel(se) > 1 && all(se(:) == se(1)), se=se(1); end
if numel(sm) > 1 && all(sm(:) == sm(1)), sm=sm(1); end

if dim > 0
  tosort=dim;
else
  tosort=1:ndims(b);
end
was_sorted=0;
myisvector = @(c)length(c) == numel(c);

for index=tosort
  x = getaxis(a, index);
  [x, sorti] = sort(x(:), index, mode);
  if ~isequal(sorti, 1:size(b, index)) && ~all(sorti == sorti(1))
    S.type = '()';
    S.subs = {};
    nx = ndims(x);
    if myisvector(x), nx=1; end
    if ndims(b) == nx
      S.subs={ sorti };
    else
      for j=1:ndims(b), 
        if j ~= index, S.subs{j}=':';
        else           S.subs{j}=sorti; end
      end
    end

    sd =subsref(sd, S);

    try
      se =subsref(se, S);
    catch
      se=[];
    end
    try
      sm =subsref(sm, S);
    catch
      sm=[];
    end
    setaxis(b, index, x);
    was_sorted=1;
  end
end
if was_sorted
  b = setalias(b, 'Signal', sd); clear sd
  label(b, 'Signal', [ 'sort(' sl ')' ]); 
  b = setalias(b, 'Error',  se);
  b = setalias(b, 'Monitor',sm);
  b.Command=cmd;
  b = history(b, mfilename, b, dim, mode);
end

