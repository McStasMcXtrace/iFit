function h = scatter3(a, option)
% h = scatter(s,option) : Plot a 2D/3D object as separate colored points
%
%   @iData/scatter function to plot a 2D or 3D object
%     2D and 3D objects are shown as separate colored points
%     The slice(a) method opens the interactive sliceomatic 3D viewer.
%   scatter(a, 'filled') produces a bubble plot which symbols are coloured and
%     proportional in size to intensity. Symbols are filled circles.
%   scatter(a, 'bubble') produces a bubble plot which symbols are coloured and
%     proportional in size to intensity. Symbols are empty circles.
%
% input:  s: object or array (iData)
%         option: global option for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel, colorbar
%                 axis tight, axis auto, view2, view3, hide_axes
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
% output: h: graphics object handles (cell)
% ex:     scatter(iData(peaks)); scatter(iData(flow));
%
% Version: $Revision$
% See also iData, iData/plot

if nargin ==1
	option='';
end
h = plot(a, [ 'scatter3 ' option ]);
if ndims(a) <= 2
  view([90 0])
end


