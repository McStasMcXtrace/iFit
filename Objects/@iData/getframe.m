function frame = getframe(a, dim, options)
% frame = getframe(s, dim, options) : create an iData object thumbnail/frame matrix
%
%   @iData/getframe function to create iData frames/thumbnails.
%   Such thumbnails may be used to create icons on GUIs
%   Note: a figure will flicker on screen.
%
% input:  s:   object or array (iData)
%         dim: dimension of the thumbnail. When specified, the frame has no labels, ticks, ...
%              to serve as a thumbnail
%         options: specific plot options, default is 'axis tight'
%
% output: frame: frame/thumbnail
% ex:     f=getframe(a); image(f.cdata);
%
% Version: $Date$
% See also iData, iData/plot, getframe, image, imwrite

persistent gl

if nargin < 2, dim=0; end
if nargin < 3, options=''; end

% test if openGL is installed correctly
if ~exist('gl') || isempty(gl)
  gl = opengl('data');
end

if numel(a) > 1
  frame = cell(size(a));
  for index=1:numel(a)
    frame{index} = getframe(a(index), dim);
  end
  return
end

if isfield(gl, 'Version') && isempty(gl.Version)
  % installed incorrectly: will lead to SEGFAULT
  frame = struct();
  return
end

f=figure('menubar','none','toolbar','none','visible','off');
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
plot(a, options); 
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
if exist('hardcopy')
  try
    frame=im2frame(hardcopy(f, '-Dzbuffer', '-r0'));
    usegetframe = 0;
  end
end 
if usegetframe
  % force figure to be 'onscreen'
  movegui(f);
  set(f,'renderer','zbuffer');
  frame=getframe(f);
end

delete(f);

