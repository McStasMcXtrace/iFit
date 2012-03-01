function h = surf(a, varargin)
% h = surf(s,option) : Plot a 2D/3D object as surface
%
%   @iData/surf function to plot a 2D or 3D object
%     2D objects are shown as a surface
%     3D objects are shown as an isosurface with median value
%     The slice(a) method opens the interactive sliceomatic 3D viewer.
%   surf(a,colormap1, b, colormap2, ..., options)
%     plots all surfaces with a specific colormap for each
%   surf(a,b, ... option)
%     same as above with a default choice of colormaps
%
% input:  s: object or array (iData)
%         option: global option for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel
%                 axis tight, axis auto, view2, view3
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
% output: h: graphics object handles (cell)
% ex:     surf(iData(peaks)); surf(iData(flow));
%
% Version: $Revision: 1.5 $
% See also iData, iData/plot

if numel(varargin) >= 1 && (isa(varargin{1},'iData') || isnumeric(varargin{1}))
  h = doublesurf(a, varargin{:}, 'surf');
elseif isempty(varargin)
  h = plot(a, [ 'surf ' ]);
else
  h = plot(a, [ 'surf ' varargin{end} ]);
end

