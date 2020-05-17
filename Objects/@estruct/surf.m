function h = surf(a, option)
% SURF   3-D colored surface.
%   H = SURFL(S) plot a 2D/3D object as a surface.
%     2D objects are shown as a surface.
%     3D objects are shown as an isosurface with median value.
%   The plot graphic object handle H is returned.
%   You may also use the SLICE(S) method to open an interactive 3D viewer.
%
%   H = SURFL(...,'OPTION') specifies plot options for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel
%                 axis tight, axis auto, view2, view3
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%
% Example: h=surfc(estruct(flow)); tf=ishandle(h); delete(gcf); tf
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot, estruct/colormap, estruct/caxis, estruct/mesh
%          estruct/slice, estruct/contour3, estruct/contourf, estruct/mesh
%          estruct/contour, estruct/slice, estruct/plot3, estruct/surfl, estruct/surfc

if nargin ==1
	option='';
end
h = plot(a, [ 'surf ' option ]);


