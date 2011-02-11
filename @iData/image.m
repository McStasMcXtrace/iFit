function h = image(r, g, b, option)
% h = image(r,g,b, option) : Plot up to 3 images overlayed on the RGB channels
%
%   @iData/image function to plot up to 3 data sets on the RBG channels of an image
%     the 3 objects should be 2D objects. Each object is re-scaled within 0 and 1
%     interval, so that each color channel is fully used. However, the 'norm'
%     keyword in options will scale the color channels within the global min and 
%     max Signal values, allowing to compare the relative intensity of the objects.
%
% input:  r: 2D object on the Red channel (iData or empty)
%         g: 2D object on the Green channel (iData or empty)
%         b: 2D object on the Blue channel (iData or empty)
%         option: global option for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel
%                 axis tight, axis auto, view2, view3, hide_axes
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%                 norm (will use the same color-scale)
% output: h: graphics object handle
% ex:     image(iData(peaks),'hide axes');
%
% Version: $Revision: 1.2 $
% See also iData, iData/plot

if nargin < 4, option = []; end
if nargin < 3, b=''; end
if nargin < 2, g=''; end

if isempty(r), r=iData; end
if isempty(g), g=iData; end
if isempty(b), b=iData; end


h = [];
a = [ r g(1) b(1) ];

if length(a(:)) > 3
  a=a(1:3);
end

% need at least one 2D object
if ~all(0 <= ndims(a) & ndims(a) <= 2), return; end
if ~any(ndims(a) == 2), return; end 

% compute larger axes
u = union(a); u=u(1);
x=getaxis(u,2);
y=getaxis(u,1);

% normalize to monitor and compute min/max
minv=Inf; maxv=-Inf;
for index=1:length(a)
  this = a(index);
  if ~isempty(this)
    this = interp(this,u);
    s = getaxis(this,0);
    if strfind(option,'norm')
      minv = min(minv, min(s(:))); 
      maxv = max(maxv, max(s(:)));
    end
    
    if     index==1, l = 'Red';
    elseif index==2, l = 'Green';
    else             l = 'Blue';
    end
    set(this, 'Monitor', 1, 'Signal', s);
    a(index) = this;
    % this now contains a valid object. We keep it for 1D to 2D extension
    if ndims(this) == 2
      set(this, 'Signal', s);
      this2D = this;
    end
  end
end

% from there, this2D has the right dimension from a previous
u=zeros(size(this2D,1),size(this2D,2),3); 

lab = '';
% normalize within [min,max]
for index=1:length(a)
  this = a(index);
  if ndims(this) == 1, this=this2D+this; end  % extend to 2D
  if ~isempty(this)
    if     index==1, l = 'Red';
    elseif index==2, l = 'Green';
    else             l = 'Blue';
    end
    lab = [ lab l ': ' this.Title ' ' ];
    s = getaxis(this,0);
    if strfind(option,'norm')
      s=s-minv; s=s/(maxv-minv);
      max(s(:))
    else
      s=s-min(s(:)); s=s/max(s(:));
    end
    u(:,:,index)=s;
  end
end

% create object to display with plot('image')
this2D.Data.cdata = u;
set(this2D, 'Signal', 'Data.cdata', 'Title', strtrim(lab));
setaxis(this2D, 2, x);
setaxis(this2D, 1, y);

% assemble the new object
h=plot(this2D, option);

