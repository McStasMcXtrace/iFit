function frame = getframe(a, dim, options)
% GETFRAME Get movie frame.
%   FRAME = GETFRAME(S) returns a movie frame/thumbnails of object S. 
%   Such thumbnails may be used to create icons on GUIs.
%   Note: a figure will flicker on screen. FRAME is a structure 
%   having the fields "cdata" and "colormap" which contain the
%   the image data in a uint8 matrix and the colormap in a double
%   matrix. F.cdata will be Height-by-Width-by-3 and F.colormap  
%   will be empty on systems that use TrueColor graphics.  
%
%   FRAME = GETFRAME(S, DIM) specifies the pixel dimensions of the image frame.
%   When specified, the frame has no labels, ticks, ...
%
%   FRAME = GETFRAME(S, DIM, 'OPTION') specific plot options.
%   Default is 'axis tight'.
%
% Example: f=getframe(estruct(peaks)); image(f.cdata); close(gcf); isstruct(f)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot, getframe, image, imwrite

if nargin < 2, dim=0; end
if nargin < 3, options=''; end

if numel(a) > 1
  frame = cell(size(a));
  for index=1:numel(a)
    frame{index} = getframe(a(index), dim);
  end
  return
end

% select a 'fast' and reliable renderer. OpenGL export often lead to black PDF/PNG...
% we use Zbuffer when available (transparently replaced by OpenGL for > R2014b)
if ndims(a) >= 2
  renderer = 'zbuffer';
else
  renderer = 'painters';
end

f=figure('menubar','none','toolbar','none','visible','off','Renderer',renderer);
% put window out of sight
p=get(f,'Position'); p(1:2) = [-1000 -1000];
set(f,'Position',p);
if dim
  % resample a with a coarse grid if needed
  need_resample=0;
  a_axes = cell(1,ndims(a));
  for index=1:ndims(a)
    s = size(a,index);  % loads object axes, or 1:end if not defined 
    if s > dim
      a_axes{index} = round(linspace(1,s,dim));
      need_resample = 1;
    else
      a_axes{index} = 1:s;
    end
  end
  % fast resampling only using indexes
  if need_resample
    S.type='()';
    S.subs=a_axes;
    a = subsref(a, S);
  end
end

% plot the data
plot(a, [ renderer ' ' options ]); 
if ndims(a) <= 2, view(2); end
drawnow
if dim
  legend off;
  xlabel(''); ylabel(''); title('');
  set(gcf,'Position',[50,50,dim,dim]);
  set(gca,'xtick',[]);
  set(gca,'ytick',[]) ;
  axis tight;
end

% extract frame
filename = [ tempname '.avi' ]; 
usegetframe = 1;
set(f,'renderer',renderer);
if exist('hardcopy')
  try
    frame=im2frame(hardcopy(f, [ '-Dzbuffer' ], '-r0'));
    usegetframe = 0;
  end
end 
if usegetframe
  % force figure to be 'onscreen'
  movegui(f);
  frame=getframe(f);
end

delete(f);

