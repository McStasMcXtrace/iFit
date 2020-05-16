function h = surfc(a, option)
% SURFC  Combination surf/contour plot.
%   H = SURFC(S) plot a 2D/3D object as surface and a contour plot
%   is drawn beneath the surface.
%   The plot graphic object handle H is returned.
%
%   H = SURFC(..., OPTION) specifies plot options for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel
%                 axis tight, axis auto, view2, view3
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%
% Example: h=surfc(estruct(flow)); tf=ishandle(h); delete(gcf); tf
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot
%          estruct/slice, estruct/contour3, estruct/contourf, estruct/mesh
%          estruct/contour, estruct/slice, estruct/plot3, estruct/surfl, estruct/surf

if nargin ==1
	option='';
end
h = plot(a, [ 'surfc ' option ]);



