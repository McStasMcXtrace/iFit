function a = scatter3(a, option)
% h = scatter3(s,option) : Plot a 2D/3D object as separate colored points
%
%   @iData/scatter3 function to plot a 2D or 3D object
%     2D and 3D objects are shown as separate colored points
%     The slice(a) method opens the interactive sliceomatic 3D viewer.
%
% input:  s: object or array (iData)
%         option: global option for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel, colorbar
%                 axis tight, axis auto, view2, view3, hide_axes
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
% output: h: graphics object handles (cell)
% ex:     scatter3(iData(peaks)); scatter3(iData(flow));
%
% Version: $Revision: 1.4 $
% See also iData, iData/plot

if nargin ==1
	option='';
end
h = plot(a, [ 'scatter3 ' option ]);



