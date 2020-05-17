function h = scatter3(a, option)
% SCATTER Scatter/bubble 2D-3D plot.
%   SCATTER(A) displays colored circles at the locations. You may also use
%   the method SLICE(A) to open the interactive sliceomatic 3D viewer.
%   SCATTER and SCATTER3 are equivalent.
%
%   H = SCATTER(A,'OPTION') specifies plot options for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel, colorbar
%                 axis tight, axis auto, view2, view3, hide_axes
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%
%   SCATTER(A, 'filled') produces a bubble plot which symbols are coloured and
%   proportional in size to intensity. Symbols are filled circles.
%
%   SCATTER(A, 'bubble') produces a bubble plot which symbols are coloured and
%   proportional in size to intensity. Symbols are empty circles.
%
% Example: h=scatter(estruct(peaks));  tf=ishandle(h); delete(gcf); tf
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot, estruct/slice, estruct/surf

if nargin ==1
	option='';
end
h = plot(a, [ 'scatter3 ' option ]);
if ndims(a) <= 2
  view([90 0])
end


