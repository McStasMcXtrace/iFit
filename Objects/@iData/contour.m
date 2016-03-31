function h = contour(a, option, varargin)
% h = contour(s,option) : Plot a 2D/3D object as contour
%
%   @iData/contour function to plot a 2D or 3D object
%     2D objects are shown as a contour
%   contour(a, option, levels) specifies the number of coutour lines to plot.
%
% input:  s: object or array (iData)
%         option: global option for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel, colorbar, hide_axes
%                 axis tight, axis auto, view2, view3
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
% output: h: graphics object handles (cell)
% ex:     contour(iData(peaks)); contour(iData(flow));
%
% Version: $Date$
% See also iData, iData/plot, iData/contour3, iData/contourf

if nargin <=1
	option='';
end

if ischar(option)
  h = plot(a, [ 'contour ' option ], varargin{:});
else
  h = plot(a, 'contour', option , varargin{:});
end



