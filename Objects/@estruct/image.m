function h = image(r, g, b, option)
% IMAGE  Display image.
%   H = IMAGE(R,G,B) plots up to 3 images overlayed as RGB channels.
%   The 3 objects should be 2D objects. Each object is re-scaled within 0 and 1
%   interval, so that each color channel is fully used. 
%
%   H = IMAGE(R,G,B, 'OPTION') specifies plot options for 2D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel
%                 axis tight, axis auto, view2, view3, hide_axes
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%                 norm (will use the same color-scale)
%
%   The 'norm' keyword in OPTION will scale the color channels within the global 
%   min and  max Signal values, allowing to compare the relative intensity of 
%   the objects.
%
% Example: image(estruct(peaks),[],[], 'hide axes');
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot

if nargin < 4, option = ''; end
if nargin < 3, b=[]; end
if nargin < 2, g=[]; end

if ischar(g), options=g; g=[]; end
if ischar(b), options=b; b=[]; end

a = [r g b ];
if length(a) == 1 && ndims(a)==3 && size(a,3)==3
  % the dataset has 3 color channels
  r=a; g=a; b=a; s=get(a,'Signal');
  set(r,'Signal',s(:,:,1));
  set(g,'Signal',s(:,:,2));
  set(b,'Signal',s(:,:,3));
  a = [ r g b ];
end

if isempty(r), r=estruct; end
if isempty(g), g=estruct; end
if isempty(b), b=estruct; end

h = [];

if numel(a) > 3
  a=a(1:3);
end

% need at least one 2D object
if ~all(0 <= ndims(a) & ndims(a) <= 2), return; end
if ~any(ndims(a) == 2), return; end 

% compute larger axes
u = union(a); u=u(1);
x=getaxis(u,2); x=double(x);
y=getaxis(u,1); y=double(y);

% normalize to monitor and compute min/max
minv=Inf; maxv=-Inf;
for index=1:numel(a)
  this = a(index);
  if ~isempty(this)
    this = interp(this,u);
    s = getaxis(this,0); s=double(s);
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
for index=1:numel(a)
  this = a(index);
  if ndims(this) == 1, this=this2D+this; end  % extend to 2D
  if ~isempty(this)
    if     index==1, l = 'Red';
    elseif index==2, l = 'Green';
    else             l = 'Blue';
    end
    lab = [ lab l ': ' this.Name ' ' ];
    s = getaxis(this,0);
    if strfind(option,'norm')
      s=s-minv; s=s/(maxv-minv);
    else
      if ~isinf(min(s(:))) && ~isinf(max(s(:)))
        s=s-min(s(:)); s=s/max(s(:));
      else
        i=find(s>=0 & isfinite(s));
        if isempty(i), i=find(isfinite(s)); end
        j=find(s<0 | ~isfinite(s));
        if ~isempty(i) && ~isempty(j), s(j)=min(s(i)); s=s/max(s(i)); end
      end
    end
    u(:,:,index)=s;
  end
end

% create object to display with plot('image')
this2D.Data.cdata = u;
setalias(this2D, 'Signal', 'Data.cdata');
this2D.Name = strtrim(lab);
setaxis(this2D, 2, x);
setaxis(this2D, 1, y);
axescheck(this2D,'force');
% assemble the new object
h=plot(this2D, [ option ' image' ]);

