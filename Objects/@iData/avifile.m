function f=avifile(a, filename, options)
% f = avifile(s, filename, options) : save iData object into a movie
%
%   @iData/avifile function to save 2D/3D data sets as AVI Movie


  if nargin < 2, filename = []; end
  if nargin < 3, options = ''; end
 
  f=saveas(a, filename, 'avi',options);
