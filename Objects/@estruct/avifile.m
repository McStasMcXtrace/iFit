function f=avifile(a, filename, options)
% AVIFILE Create a new AVI file by rotating the object 2D/3D view.
%   F = AVIFILE(S) creates a movie file named from the object Tag.
%
%   F = AVIFILE(S, FILENAME) saves object into given movie file.
%
%   F = AVIFILE(S, FILENAME, 'OPTIONS') specifies plot options, such as:
%     tight, hide_axes, interp, view2, view3, transparent, light, clabel, colorbar, 
%     whole. Default is 'view3 axis tight'.
%
% Example: f=avifile(estruct(peaks)); tf=~isempty(dir(f)); delete(f); tf
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/saveas

  if nargin < 2, filename = []; end
  if nargin < 3, options = ''; end
 
  f=saveas(a, filename, 'avi',options);
