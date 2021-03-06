function h = plot3(a, option)
% PLOT3   Plot a 2D/3D object as separate points or semi-transparent volume.
%   PLOT3(A) plots lines and points in 3-D space.
%     2D objects are shown as separate points in one single color
%     3D objects are shown as a semi-transparent volume
%   The SLICE(A) method opens the interactive sliceomatic 3D viewer.
%
%   H = PLOT3(A,'OPTION') specifies plot options for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel
%                 axis tight, axis auto, view2, view3
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%
% Example: h=plot3(iData(flow)); tf=ishandle(h); delete(gcf); tf
% Version: $Date$ $Version$ $Author$
% See also iData, iData/plot, iData/slice, iData/contour3, iData/contourf
%          iData/contour, iData/surf, iData/slice, iData/surfl, iData/surfc, iData/mesh

if nargin ==1
	option='';
end
h = plot(a, [ 'plot3 ' option ]);



