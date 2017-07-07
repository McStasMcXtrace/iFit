function f=imwrite(a, filename, format)
% f = imwrite(s, filename, format) : save iData object into an image
%
%   @iData/saveas function to a data set as an image (PNG is default)

  if nargin < 2, filename = []; end
  if nargin < 3, format = 'png'; end
  
  f = saveas(a, filename, format);
