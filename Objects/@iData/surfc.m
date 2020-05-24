function h = surfc(a, option)
% SURFC  Combination surf/contour plot.
%   H = SURFC(S) plot a 2D/3D object as surface and a contour plot
%   is drawn beneath the surface.
%     2D objects are shown as a surface.
%     3D objects are shown as an isosurface with median value.
%   The plot graphic object handle H is returned.
%
%   H = SURFC(..., 'OPTION') specifies plot options for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel
%                 axis tight, axis auto, view2, view3
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%
% Example: h=surfc(iData(flow)); tf=ishandle(h); delete(gcf); tf
% Version: $Date$ $Version$ $Author$
% See also iData, iData/plot
%          iData/slice, iData/contour3, iData/contourf, iData/mesh
%          iData/contour, iData/slice, iData/plot3, iData/surfl, iData/surf

if nargin ==1
	option='';
end
h = plot(a, [ 'surfc ' option ]);



