function h = colormap(varargin)
% colormap: Produce surfaces from iData 2D objects with different colormaps
%
% colormap(z1,cm1, z2, cm2, ..., option)
%   Produces surfaces with different colormaps.
%   iData 2D object z1 is plotted as surface 1, which is coloured with map cm1 (Nx3)
%   ...
%   options specify the type of plot rendering.
%
% colormap(z1, z2, ...)
% Same as above with a set of default colormaps
%
% Any string argument will be interpreted as an 'option' for the 2D iData/plot
%                 surf, mesh, contour, contour3, surfc, surfl, contourf
%                 plot3, scatter3 (colored points), stem3, pcolor, waterfall
%                 flat, interp, faceted (for shading), view2, view3
%                 transparent, light, clabel, colorbar, shifted (overlayed 2D)
%                 axis tight, axis auto, hide_axes (compact layout)
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%                 whole (do not reduce large object size for plotting)
%
% Available colormaps are:
%    hsv        - Hue-saturation-value color map.
%    hot        - Black-red-yellow-white color map.
%    gray       - Linear gray-scale color map.
%    bone       - Gray-scale with tinge of blue color map.
%    copper     - Linear copper-tone color map.
%    pink       - Pastel shades of pink color map.
%    white      - All white color map.
%    flag       - Alternating red, white, blue, and black color map.
%    lines      - Color map with the line colors.
%    colorcube  - Enhanced color-cube color map.
%    vga        - Windows colormap for 16 colors.
%    jet        - Variant of HSV.
%    prism      - Prism color map.
%    cool       - Shades of cyan and magenta color map.
%    autumn     - Shades of red and yellow color map.
%    spring     - Shades of magenta and yellow color map.
%    winter     - Shades of blue and green color map.
%    summer     - Shades of green and yellow color map.

% Copyright (c) 2001 by OPTI-NUM solutions
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/29419 ; 28/02/2012

z=[]; cm = []; options='';
for index=1:length(varargin)
  this = varargin{index};
  if     isa(this, 'iData') && all(ndims(this)==2), z = [ z this ]; 
  elseif isnumeric(this) && size(this,2)==3
    cm = [ cm {this} ];
  elseif ischar(this), options=[ options ' ' this ]; 
  end
end
cm= [ cm {hot} {cool} {winter} {spring} {jet} {hsv} {bone} {copper} {pink} ]; % default colormaps at the end in case too few are defined

% Build the actual colormap by catenation
cmap = cat(1, cm{:});

% Now we make up the color indices.
sumcm = 0;
ci= {};
for index=1:numel(z)
  % compute local colormap so that it matches the object values
  this = double(z(index)); 
  if any(this < 0) this = this-min(this(:)); end
  zfloor=linspace(min(this(:)),max(this(:)), size(cm{index},1));
  cind  =zeros(size(this));
  % count elements which are lower than each colormap value zfloor
  % sum up all colormaps on top of each other
  for k=1:prod(size(this)), 
    cind(k)=sum(zfloor<=this(k))+sumcm; 
  end
  ci{index} = cind;
  sumcm = sumcm+size(cm{index},1);
end

% And now we make the surfaces:
h = plot(z, options);
for index=1:numel(z)
  this = z(index);
  try
    set(h(index), 'CDataMapping','direct', 'CData', ci{index});
  end
  if index==1
    title(title(this));
    xlabel(xlabel(this));
    if ndims(this) >= 2, ylabel(ylabel(this)); end
  end
end

colormap(cmap);

