function f=imwrite(a, filename, format)
% f = imwrite(s, filename, format) : save iFunc object into an image
%
%   @iFunc/imwrite function to save mode as an image (PNG is default)

  if nargin < 2, filename = []; end
  if nargin < 3, format = 'png'; end
  
  f = saveas(a, filename, format);
