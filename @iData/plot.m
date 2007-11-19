function h=plot(a, method)
% h = plot(s, method) : plot iData object
%
%   @iData/plot function to plot data sets
%   This function plot the signal of the object as a function of the defined axes.
%   The plot is an errorbar plot for 1D data y=f(x), a surface mesh for 2D z=f(x,y),
%     and an isosurface for 3D c=f(x,y,z) data. Data specified as series of 
%     coordinate points are plotted using a plot3-type rendering. Further
%     dimensionalities are not handled.
%
% input:  s: object or array (iData)
%         method: optional type of plot to render for 2D and 3D views, within
%                 surf, mesh, contour, contour3, surfc, surfl, contourf, stem3
%                 flat, interp, faceted, transparent, light, clabel
%                 plot3, scatter3
% output: h: graphics object handles (cell)
% ex:     plot(iData(rand(10), 'surfc interp transparent');
%
% See also iData, interp1, interpn, ndgrid, iData/setaxis, iData/getaxis
% Uses fscatter3 from Felix Morsdorf, Jan 2003, Remote Sensing Laboratory Zuerich

% private function:
%   fscatter3: Felix Morsdorf, Jan 2003, Remote Sensing Laboratory Zuerich
%   vol3d:     Joe Conti, 2004
if nargin == 1, method=''; end
if length(a) > 1
  h = cell(size(a));
  for index=1:length(a(:))
    h{index} = plot(a(index), method);
  end
  h = reshape(h, size(a));
  return
end

switch ndims(a)
case 0
  h=[];
case 1  % vector type data (1 axis + signal) -> plot
  [x, xlab] = getaxis(a,1);
  [y, ylab] = getaxis(a,0);
  e         = get(a,'Error');
  m         = get(a,'Monitor');
  if not(all(m == 1) | all(m == 0)) & a.PerMonitor, 
    y = y./m; e=e./m; ylab = [ylab ' per monitor' ];
  end
  h = errorbar(x,y,e);
case 2  % surface type data (2 axes+signal) -> surf or plot3
  [x, xlab] = getaxis(a,1);
  [y, ylab] = getaxis(a,2);
  [z, zlab] = getaxis(a,0);
  m         = get(a,'Monitor');
  if not(all(m == 1) | all(m == 0)) & a.PerMonitor, 
    z = z./m; zlab = [zlab ' per monitor' ];
  end
  if isvector(a) == 2 % plot3
    if (strfind(method,'plot3'))
      h = plot3(x,y,z);
    else
      h=fscatter3(x,y,z,z);
    end
  else                % surf
    C = [];
    if (strfind(method,'contour3'))
      [C,h] =contour3(x,y,z);
    elseif (strfind(method,'contourf'))
      [C,h]=contourf(x,y,z);
    elseif (strfind(method,'contour'))
      [C,h]=contour(x,y,z);
    elseif (strfind(method,'surfc'))
      h=surfc(x,y,z);
    elseif (strfind(method,'surfl'))
      h=surfl(x,y,z);
    elseif (strfind(method,'mesh'))
      h=mesh(x,y,z);
    elseif (strfind(method,'stem3'))
      h=stem3(x,y,z);
    else
      h=surf(x,y,z);
    end

    if ~isempty(C) & strfind(method,'clabel')
      clabel(C,h);
    end
  end
  zlabel(zlab);
case 3
  [x, xlab] = getaxis(a,1);
  [y, ylab] = getaxis(a,2);
  [z, zlab] = getaxis(a,3);
  [c, clab] = getaxis(a,0);
  m         = get(a,'Monitor');
  if not(all(m == 1) | all(m == 0)) & a.PerMonitor, c = c./m; clab = [clab ' per monitor' ]; end
  if isvector(a) == 3 % plot3-like
    h=fscatter3(x,y,z,c);
  else
    if strfind(method, 'surf')
      h = vol3d('cdata',c,'texture','3D','xdata',x,'ydata',y,'zdata',z);
      alphamap('vdown');
      vol3d(h);
      h = h.handles;
    else
      a = interp(a,'grid');
      isosurface(x,y,z,c,median(c(:)));
      h = findobj(gca,'type','patch');
      hold off
    end
  end;
  zlabel(zlab);
otherwise
  iData_private_warning(mfilename, [ 'plotting of ' num2str(ndims(a)) '-th dimensional data is not implemented from ' a.Tag ]);
end

if (strfind(method,'flat'))
  shading flat
elseif (strfind(method,'interp'))
  shading interp
elseif (strfind(method,'faceted'))
  shading faceted
end
if (strfind(method,'transparent') | strfind(method,'alpha'))
  alpha(0.7);
end
if (strfind(method,'light'))
  light;
end

% add a UIcontextMenu so that right-click gives info about the iData plot
% also make up title string
uicm = uicontextmenu;
uimenu(uicm, 'Label', [ 'Data ' a.Tag ]);
uimenu(uicm, 'Label', [ num2str(ndims(a)) 'D object ' mat2str(size(a)) ]);
T=a.Title; if iscell(T), T=T{1}; end
T=deblank(T);
if length(T) > 23, T=[ T(1:20) '...' ]; end
titl=T;
uimenu(uicm, 'Label', [ 'Title: "' T '"' ]);
T=a.Source;
if length(T) > 23, T=[ '...' T(end-20:end) ]; end
cmd = a.Command{end};
if length(cmd) > 23, cmd = [ cmd(1:20) '...' ]; end
titl =[ titl ' <' T '> (' a.Tag ':' cmd ')' ];
uimenu(uicm, 'Label', [ 'Source: <' T '>' ]);
uimenu(uicm, 'Label', [ 'Cmd: ' cmd ]);
uimenu(uicm, 'Separator','on','Label','Reset View', 'Callback','view(0,90);lighting none;alpha(1);shading faceted;axis auto');
uimenu(uicm, 'Label','Toggle grid', 'Callback','grid');
if ndims(a) >= 2
  uimenu(uicm, 'Label','Smooth View','Callback', 'shading interp;');
  uimenu(uicm, 'Label','Add Light','Callback', 'light;lighting phong;');
  uimenu(uicm, 'Label','Transparency','Callback', 'alpha(0.7);');
end
% attach menu to plot
set(h, 'UIContextMenu', uicm);
set(h, 'Tag', char(a));
set(gcf, 'Name', char(a));

% labels
xlabel(xlab);
ylabel(ylab);
if ndims(a) == 3
  title({ clab ,titl });
else
  title(titl);
end

% ============================================================================
function [h] = fscatter3(X,Y,Z,C,cmap);
% [h] = fscatter3(x,y,z,C,cmap);
% Plots point cloud data in cmap color classes and 3 Dimensions,
% much faster and very little memory usage compared to scatter3 !
% x,y,z,C are vectors of the same length
% with C being used as index into colormap (can be any values though)
% cmap is optional colourmap to be used
% h are handles to the line objects

% Felix Morsdorf, Jan 2003, Remote Sensing Laboratory Zuerich

if nargin == 3
  C = Z;
end
if nargin <= 4
  numclass = 256; % Number of color classes
  cmap = hsv(256);
end
if nargin == 5
  numclass = max(size(cmap));
  if numclass == 1
    cmap = hsv(256);
    numclass = 256;
  end  
end

% avoid too many calculations

mins = min(C);
maxs = max(C);
minz = min(Z);
maxz = max(Z);
minx = min(X);
maxx = max(X);
miny = min(Y);
maxy = max(Y);

% construct colormap :

col = cmap;

% determine index into colormap

ii = round(interp1([floor(mins) ceil(maxs)],[1 numclass],C));
hold on
colormap(cmap);

% plot each color class in a loop

marker = max(4, 7-log10(length(X(:))));

k = 0;
for j = 1:numclass
  jj = find(ii == j);
  if ~isempty(jj)
    k = k + 1;
    h(k) = plot3(X(jj),Y(jj),Z(jj),'.','color',col(j,:), ...
		 'markersize',marker);
  end  
end
hold off

% ============================================================================
function [model] = vol3d(varargin)
%H = VOL3D Volume render 3-D data. 
% VOL3D uses the orthogonal plane 2-D texture mapping technique for 
% volume rending 3-D data in OpenGL. Use the 'texture' option to fine 
% tune the texture mapping technique. This function is best used with
% fast OpenGL hardware.
%
% H = vol3d('CData',data) Create volume render object from input 
%                         3-D data. Use interp3 on data to increase volume
%                         rendering resolution. Returns a struct 
%                         encapsulating the pseudo-volume rendering object. 
%
% vol3d(...'Parent',axH) Specify parent axes
%
% vol3d(...,'texture','2D')  Default. Only render texture planes parallel
%                            to nearest orthogonal viewing plane. Requires
%                            doing vol3d(h) to refresh if the view is
%                            rotated (i.e. using cameratoolbar).
%
% vol3d(...,'texture','3D')  Render x,y,z texture planes simultaneously. 
%                            This avoids the need to refresh the view but 
%                            requires faster OpenGL hardware peformance.
%
% vol3d(H)  Refresh view. Updates rendering of texture planes 
%           to reduce visual aliasing when using the 'texture'='2D'
%           option.
%
% NOTES
% Use vol3dtool for editing the colormap and alphamap. 
% Adjusting these maps will allow you to explore your 3-D volume 
% data at various intensity levels. See documentation on 
% alphamap and colormap for more information.
%
% Use interp3 on input date to increase/decrease resolution of data
%
% Examples:
%
% % Visualizing fluid flow
% v = flow(50); [x,y,z,v]=flow;
% h = vol3d('cdata',v,'texture','2D');
% h = vol3d('cdata',v,'texture','2D','xdata',x,'ydata',y,'zdata',z);
% view(3); 
% % Update view since 'texture' = '2D'
% vol3d(h);  
% alphamap('rampdown'), alphamap('decrease'), alphamap('decrease')
% 
% % Visualizing MRI data
% load mri.mat
% D = squeeze(D);
% h = vol3d('cdata',D,'texture','2D');
% view(3); 
% % Update view since 'texture' = '2D'
% vol3d(h);  
% axis tight;  daspect([1 1 .4])
% alphamap('rampup');
% alphamap(.06 .* alphamap);
%
% See also vol3dtool, alphamap, colormap, opengl, isosurface

% Copyright Joe Conti, 2004

if isstruct(varargin{1})
    model = varargin{1};
    if length(varargin) > 1
       varargin = {varargin{2:end}};
    end
else
    model = localGetDefaultModel;
end


if length(varargin)>1
  for n = 1:2:length(varargin)
    switch(lower(varargin{n}))
        case 'cdata'
            model.cdata = varargin{n+1};
        case 'parent'
            model.parent = varargin{n+1};
        case 'texture'
            model.texture = varargin{n+1};
        case 'xdata'
            model.xdata = varargin{n+1};
            if length(model.xdata) > 2
              model.xdata = [ model.xdata(1) model.xdata(end) ];
            end
        case 'ydata'
            model.ydata = varargin{n+1};
            if length(model.ydata) > 2
              model.ydata = [ model.ydata(1) model.ydata(end) ];
            end
        case 'zdata'
            model.zdata = varargin{n+1};
            if length(model.zdata) > 2
              model.zdata = [ model.zdata(1) model.zdata(end) ];
            end
    end
    
  end
end

if isempty(model.parent)
    model.parent = gca;
end

% choose default 3-D view
ax = model.parent;
axis(ax,'vis3d');
axis(ax,'tight');

[model] = local_draw(model);


%------------------------------------------%
function [model] = localGetDefaultModel

model.cdata = [];
model.xdata = [];
model.ydata = [];
model.zdata = [];
model.parent = [];
model.handles = [];
model.texture = '2D';

%------------------------------------------%
function [model,ax] = local_draw(model)

cdata = model.cdata; 
siz = size(cdata);

% Define [x,y,z]data
if isempty(model.xdata)
    model.xdata = [0 siz(2)];
end
if isempty(model.ydata)
    model.ydata = [0 siz(1)];
end
if isempty(model.zdata)
    model.zdata = [0 siz(3)];
end

try,
   delete(model.handles);
end

ax = model.parent;
cam_dir = camtarget(ax) - campos(ax);
[m,ind] = max(abs(cam_dir));

is3DTexture = strcmpi(model.texture,'3D');
handle_ind = 1;

% Create z-slice
if(ind==3 | is3DTexture )    
  x = [model.xdata(1), model.xdata(2); model.xdata(1), model.xdata(2)];
  y = [model.ydata(1), model.ydata(1); model.ydata(2), model.ydata(2)];
  z = [model.zdata(1), model.zdata(1); model.zdata(1), model.zdata(1)];
  diff = model.zdata(2)-model.zdata(1);
  delta = diff/size(cdata,3);
  for n = 1:size(cdata,3)

   slice = double(squeeze(cdata(:,:,n)));
   h(handle_ind) = surface(x,y,z,'Parent',ax);
   set(h(handle_ind),'cdatamapping','scaled','facecolor','texture','cdata',slice,...
	 'edgealpha',0,'alphadata',double(slice),'facealpha','texturemap','tag','vol3d');
   z = z + delta;
   handle_ind = handle_ind + 1;
  end

end

% Create x-slice
if (ind==1 | is3DTexture ) 
  x = [model.xdata(1), model.xdata(1); model.xdata(1), model.xdata(1)];
  y = [model.ydata(1), model.ydata(1); model.ydata(2), model.ydata(2)];
  z = [model.zdata(1), model.zdata(2); model.zdata(1), model.zdata(2)];
  diff = model.xdata(2)-model.xdata(1);
  delta = diff/size(cdata,2);
  for n = 1:size(cdata,2)

   slice = double(squeeze(cdata(:,n,:)));
   h(handle_ind) = surface(x,y,z,'Parent',ax);
   set(h(handle_ind),'cdatamapping','scaled','facecolor','texture','cdata',slice,...
	 'edgealpha',0,'alphadata',double(slice),'facealpha','texturemap','tag','vol3d');
   x = x + delta;
   handle_ind = handle_ind + 1;
  end
end

  
% Create y-slice
if (ind==2 | is3DTexture)
  x = [model.xdata(1), model.xdata(1); model.xdata(2), model.xdata(2)];
  y = [model.ydata(1), model.ydata(1); model.ydata(1), model.ydata(1)];
  z = [model.zdata(1), model.zdata(2); model.zdata(1), model.zdata(2)];
  diff = model.ydata(2)-model.ydata(1);
  delta = diff/size(cdata,1);
  for n = 1:size(cdata,1)

   slice = double(squeeze(cdata(n,:,:)));
   h(handle_ind) = surface(x,y,z,'Parent',ax);
   set(h(handle_ind),'cdatamapping','scaled','facecolor','texture','cdata',slice,...
	 'edgealpha',0,'alphadata',double(slice),'facealpha','texturemap','tag','vol3d');
   y = y + delta;
   handle_ind = handle_ind + 1;
  end
end

model.handles = h;

