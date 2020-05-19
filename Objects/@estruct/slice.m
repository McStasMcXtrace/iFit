function h = slice(a, method)
% SLICE  Volumetric interactive slice plot.
%   SLICE(S) plots a 2D/3D object with slice and volume rendering.
%   The plot is obtained with Matlab Central sliceomatic. Large objects are
%   reduced to 10^6 elements. Use SLICE(S, 'whole') to use full data.
%
% Example: h=slice(estruct(flow)); tf=ishandle(h); delete(gcf); tf
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot, sliceomatic

if ndims(a) == 2
  
  fig = figure; h = plot(a);
  % we move the cuts on the right and bottom sides
  ax = get(cuthandles(1),'Parent'); fx = get(ax,'Parent');
  ay = get(cuthandles(2),'Parent'); fy = get(ay,'Parent');
  p = get(fig, 'Position');
  set(fx, 'Position', [ p(1)+p(3) p(2) p(3:4)/2 ]);
  set(fy, 'Position', [ p(1) p(2)-p(4)/2 p(3:4)/2 ]);
  figure(fig);
  return
elseif ndims(a) < 3
  warning([ mfilename ': Slice-o-matic is only available for 3D objects, but ' a.Tag ' ' a.Name '" is ' num2str(ndims(a)) '-th dimensions. Using plot instead.' ]);
  h = plot(a);
  return
elseif  isvector(a)
  warning([ mfilename ': Creating an histogram data set from data ' a.Tag ' "' a.Name '" a=hist(a);' ]);
  a = hist(a);
end
if nargin < 2
  method = '';
end

if ndims(a) > 3
  % reduce dimensions
  sz = size(a);  sz(4:end) = 1;
  warning([ mfilename ': Reducing ' num2str(ndims(a)) '-th dimensional data ' a.Tag ' "' a.Name '" to 3D with a=resize(a, ' mat2str(sz) ')' ]);
  a = resize(a, sz);
end

if prod(size(a)) > 1e6
  if isempty([ strfind(method,'whole') strfind(method,'full') ])
    warning([ mfilename ': Object ' a.Tag ' is large (numel=' num2str(prod(size(a))) ...
      ').\n\tNow rebinning for display purposes with e.g. a=reducevolume(a);' ...
      '\n\tUse e.g slice(a, ''whole'') to plot the whole data set and be able to zoom tiny regions.' ]);
    a=reducevolume(a);
  end
end

if exist('sliceomatic')
  x=unique(getaxis(a,2));
  y=unique(getaxis(a,1));
  z=unique(getaxis(a,3));
  c=getaxis(a,0);
  sliceomatic(c, x,y,z);
  h = findobj('name','Sliceomatic');
  set(h, 'Name', [ 'Sliceomatic: ' char(a) ]); 
  title(a.Name,'interpreter','none');
  xlabel(xlabel(a)); ylabel(ylabel(a)); zlabel(zlabel(a));
else
  h = plot(a);
end

