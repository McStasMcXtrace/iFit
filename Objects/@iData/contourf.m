function h = contourf(a, option, varargin)
% h = contourf(s,option) : Plot a 2D/3D object as filled contour
%
%   @iData/contourf function to plot a 2D or 3D object
%     2D objects are shown as a filled contour
%   contourf(a, option, levels) specifies the number of coutour lines to plot.
%
% input:  s: object or array (iData)
%         option: global option for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel, colorbar, hide_axes
%                 axis tight, axis auto, view2, view3
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
% output: h: graphics object handles (cell)
% ex:     contourf(iData(peaks)); contourf(iData(flow));
%
% Version: $Date$
% See also iData, iData/plot, iData/contour3, iData/contour

if nargin <=1
	option='';
end
if ischar(option)
  h = plot(a, [ 'contourf ' option ], varargin{:});
else
  h = plot(a, 'contourf', option , varargin{:});
end



