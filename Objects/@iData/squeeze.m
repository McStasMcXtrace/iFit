function a = squeeze(a)
% c = squeeze(a) : remove singleton dimensions from iData objects/arrays
%
%   @iData/squeeze returns an object B with the same elements as
%    A but with all the singleton dimensions removed.  A singleton
%    is a dimension such that size(A,dim)==1.  2-D arrays are
%    unaffected by squeeze so that row vectors remain rows.
%
% input:  a: object or array (iData)
% output: c: object or array (iData)
% ex:     c=squeeze(zeros(iData,[2 1 3]));
%
% Version: $Date$
% See also iData, iData/size, iData/resize, iData/permute, iData/reshape

d = size(a);
x = find(d <= 2);
d = d(find(d > 1));
d = [d ones(1,2-length(d))]; % Make sure siz is at least 2-D

if numel(a) == 1
  % this is a single iData object
  if ndims(a) <= 2 && length(size(a)) ==2, return; end
  s = get(a,'Signal');
  sq= squeeze(s);
  [dummy, sl] = getaxis(a, '0');  % signal definition/label
  a = set(a,'Signal',sq, [  'squeeze(' sl ')' ]);
  % check if we could update object
  if ~isequal(subsref(a,struct('type','.','subs','Signal')), sq)
    a = setalias(a, 'Signal', sq, [  'squeeze(' sl ')' ]);
  end
  for index=1:length(x)
    a=rmaxis(a, x(index));
  end
else
  % this is an iData array
  a = reshape(a, d);
  for index=1:numel(a)
      a(index) = squeeze(a(index));
  end
end

