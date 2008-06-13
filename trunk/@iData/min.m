function m = min(a,b)
% m = min(a,b) : computes the minimum value of iData object(s)
%
%   @iData/min function to compute the minimum value of data sets.
%     min(iData) returns a single value as the minimum value of the iData signal
%     min(a,b)   returns an object which signal is the lowest of a and b.
%
% input:  a: object or array (iData)
%         b: object or array (iData/double)
% output: m: minimum value (double/iData)
% ex:     b=min(a); or min(a,1)
%
% Version: $Revision: 1.3 $
% See also iData, iData/min, iData/max

if nargin == 1
  b = [];
end

% handle input iData arrays
if length(a) > 1 & isa(a,'iData')
  m = [];
  for index=1:length(a(:))
    m(index) = min(a(index), b);
  end
  m = reshape(m, size(a));
  return
end

if ~isa(a, 'iData')
  m = min(b, a);
  return
end

if isempty(b)
  m = get(a, 'Signal');
  m = min(m(:));
  return
end

% find intersection between iData objects
if isa(b, 'iData')
  [a,b] = intersect(a,b);
  m = copyobj(a);
  set(m, 'Signal', min(get(a,'Signal'), get(b,'Signal')));
  return
else
% handle iData and scalar/vector/matrix min/max
  m = copyobj(a);
  set(m, 'Signal', min(get(a,'Signal'), b));
end


