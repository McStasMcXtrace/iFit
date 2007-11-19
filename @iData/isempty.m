function y = isempty(s)
% isempty(s) : true for empty iData object
%
%   @iData/isempty true for empty iData object
%
% input:  s: object or array (iData)
% output: false(0) or true(1) whether Signal is empty in the objects
% ex :    isempty(iData)
%
% See also iData, iData/disp, iData/get, iData/size

% EF 23/09/07 iData implementation

y = zeros(size(s));

for index = 1:length(s)
  if any(size(s(index))) == 0, empty = 1;
  else                         empty = 0; end
  y(index) = empty;
end
y=uint8(y);


