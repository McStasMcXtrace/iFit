function h = contourf(a, option, varargin)
% CONTOURF Filled contour plot.
%   CONTOURF(Z) is a contour plot of 2D/3D object Z treating the values in Z
%   as heights above a plane.  A contour plot are the level curves
%   of Z for some values V. The areas between contours are filled with colors 
%   according to the Z-value for each level.  Contour regions with data values 
%   at or above a given level are filled with the color that maps to the 
%   interval. NaN's in the Z-data leave white holes with black borders in the
%   contour plot.
%
%   H = CONTOURF(..., 'OPTION') specifies plot options for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel, colorbar, hide_axes
%                 axis tight, axis auto, view2, view3
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%
%   CONTOURF(A, 'OPTION', LEVELS) specifies the number of coutour lines to plot.
%
% Example: h=contourf(iData(peaks)); tf=ishandle(h); delete(gcf); tf
% Version: $Date$ $Version$ $Author$
% See also iData, iData/plot, iData/contour3, iData/contour

if nargin <=1
	option='';
end
if ischar(option)
  h = plot(a, [ 'contourf ' option ], varargin{:});
else
  h = plot(a, 'contourf', option , varargin{:});
end



