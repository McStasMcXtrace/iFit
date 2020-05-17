function a = mesh(a, option)
% MESH   3-D mesh surface.
%   H = MESH(S) plots a 2D object as a coloured wireframe mesh.
%
%   H = MESH(...,'OPTION') specifies plot options for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel, colorbar, hide_axes
%                 axis tight, axis auto, view2, view3
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%
% Example: h=mesh(estruct(peaks)); tf=ishandle(h); delete(gcf); tf
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot
%          estruct/slice, estruct/contour3, estruct/contourf
%          estruct/contour, estruct/slice, estruct/plot3, estruct/surfl, estruct/surfc

if nargin ==1
	option='';
end
h = plot(a, [ 'mesh ' option ]);



