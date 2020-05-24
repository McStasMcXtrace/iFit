function h = contour(a, option, varargin)
% CONTOUR Contour plot.
%   CONTOUR(Z) is a contour plot of 2D/3D object Z treating the values in Z
%   as heights above a plane.  A contour plot are the level curves
%   of Z for some values V.  
%
%   H = CONTOUR(..., 'OPTION') specifies plot options for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel, colorbar, hide_axes
%                 axis tight, axis auto, view2, view3
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%
%   CONTOUR(A, 'OPTION', LEVELS) specifies the number of coutour lines to plot.
%
% Example: h=contour(iData(peaks)); tf=ishandle(h); delete(gcf); tf
% Version: $Date$ $Version$ $Author$
% See also iData, iData/plot, iData/contour3, iData/contourf

if nargin <=1
	option='';
end

if ischar(option)
  h = plot(a, [ 'contour ' option ], varargin{:});
else
  h = plot(a, 'contour', option , varargin{:});
end



