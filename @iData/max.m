function m = max(a,b)
% m = max(s) : computes the maximum value of iData object(s)
%
%   @iData/max function to compute the maximum value of data sets.
%     max(iData) returns a single value as the maximum value of the iData signal
%     max(a,b)   returns an object which signal is the lowest of a and b.
%
% input:  s: object or array (iData)
% output: m: maximum value (double/iData)
% ex:     b=max(a);
%
% See also iData, iData/max, iData/min

if nargin == 1
  b = [];
end

% handle input iData arrays
if length(a) > 1 & isa(a,'iData')
  m = [];
  for index=1:length(a(:))
    m(index) = max(a(index), b);
  end
  m = reshape(m, size(a));
  return
end

if ~isa(a, 'iData')
  m = max(b, a);
  return
end

if isempty(b)
  m = get(a, 'Signal');
  m = max(m(:));
  return
end

% find intersection between iData objects
if isa(b, 'iData')
  [a,b] = intersect(a,b);
  m = copyobj(a);
  set(m, 'Signal', max(get(a,'Signal'), get(b,'Signal')));
  return
else
% handle iData and scalar/vector/matrix max/max
  m = copyobj(a);
  set(m, 'Signal', max(get(a,'Signal'), b));
end


