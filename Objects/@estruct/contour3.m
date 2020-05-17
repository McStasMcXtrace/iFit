function h = contour3(a, option, varargin)
% CONTOUR3 3D contour plot.
%   CONTOUR3(Z) is a contour plot of 2D/3D object Z treating the values in Z
%   as heights above a plane.  A contour plot are the level curves
%   of Z for some values V. the contours are drawn at their corresponding Z level.
%
%   H = CONTOUR3(..., 'OPTION') specifies plot options for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel, colorbar, hide_axes
%                 axis tight, axis auto, view2, view3
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%
%   CONTOUR3(A, 'OPTION', LEVELS) specifies the number of coutour lines to plot.
%
% Example: h=contour3(estruct(flow)); tf=ishandle(h); delete(gcf); tf
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot, estruct/contour3, estruct/contourf

if nargin <=1
	option='';
end
if ischar(option)
  h = plot(a, [ 'contour3 ' option ], varargin{:});
else
  h = plot(a, 'contour3', option , varargin{:});
end
