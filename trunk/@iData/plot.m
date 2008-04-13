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
% Contributed code (Matlab Central): 
%   fscatter3: Felix Morsdorf, Jan 2003, Remote Sensing Laboratory Zuerich
%   vol3d:     Joe Conti, 2004

% private functions:
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

switch ndims(a) % handle different plotting methods depending on the iData dimensionality
case 0
  h=[]; return;
case 1  % vector type data (1 axis + signal) -> plot
  [x, xlab] = getaxis(a,1);
  [y, ylab] = getaxis(a,0);
  e         = get(a,'Error');
  m         = get(a,'Monitor');
  if not(all(m == 1) | all(m == 0)),
    y = y./m; e=e./m; ylab = [ylab ' per monitor' ];
  end
  if isempty(method), method='b-'; end
  if all(e == 0)
    h = plot(x,y, method);
  else
    h = errorbar(x,y,e,method);
  end
case 2  % surface type data (2 axes+signal) -> surf or plot3
  [x, xlab] = getaxis(a,1);
  [y, ylab] = getaxis(a,2);
  [z, zlab] = getaxis(a,0);
  m         = get(a,'Monitor');
  if not(all(m == 1) | all(m == 0)),
    z = z./m; zlab = [zlab ' per monitor' ];
  end
  if isvector(a) == 2 % plot3/fscatter3
    if (strfind(method,'plot3'))
      h = plot3(x,y,z);
    else
      h=fscatter3(x,y,z,z);
    end
  else                % surf and similar stuff
    C = [];
    if isvector(x) & isvector(y),
      z = z';
    end
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
case 3  % #d data sets: volumes
  % first test if this is an image
  if isfield(a.Data,'cdata')
    h=image(a.Data.cdata);
    xlab=''; ylab=''; clab='';
  else
    [x, xlab] = getaxis(a,1);
    [y, ylab] = getaxis(a,2);
    [z, zlab] = getaxis(a,3);
    [c, clab] = getaxis(a,0);
    m         = get(a,'Monitor');
    if not(all(m == 1) | all(m == 0)), c = c./m; clab = [clab ' per monitor' ]; end
    if isvector(a) == 3 | ~isempty(strfind(method, 'plot3')) % plot3-like
      h=fscatter3(x(:),y(:),z(:),c(:));
    else
      if ~isempty(strfind(method, 'surf')) | ~isempty(strfind(method, 'vol3d'))
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
  end
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
uimenu(uicm, 'Label', [ 'Data ' a.Tag ': ' num2str(ndims(a)) 'D object ' mat2str(size(a)) ]);
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
uimenu(uicm, 'Label', [ 'User: ' a.User ]);
uimenu(uicm, 'Separator','on','Label','Toggle grid', 'Callback','grid');
if ndims(a) >= 2
  uimenu(uicm, 'Label','Reset View', 'Callback','view(2);lighting none;alpha(1);shading faceted;axis tight;set(gco,''Edgecolor'',''none'');');
  uimenu(uicm, 'Label','Smooth View','Callback', 'shading interp;');
  uimenu(uicm, 'Label','Add Light','Callback', 'light;lighting phong;');
  uimenu(uicm, 'Label','Transparency','Callback', 'alpha(0.7);');
  uimenu(uicm, 'Label','Linear/Log scale','Callback', 'if strcmp(get(gca,''zscale''),''linear'')  set(gca,''zscale'',''log''); else set(gca,''zscale'',''linear''); end');
else
  uimenu(uicm, 'Label','Reset View', 'Callback','view(2);lighting none;alpha(1);shading faceted;axis tight');
  uimenu(uicm, 'Label','Linear/Log scale','Callback', 'if strcmp(get(gca,''yscale''),''linear'')  set(gca,''yscale'',''log''); else set(gca,''yscale'',''linear''); end');
end
% attach contexual menu to plot
set(h, 'UIContextMenu', uicm);
set(h, 'Tag', char(a));
set(gcf, 'Name', char(a));

% labels
if ~isempty(xlab), xlabel(xlab,'interpreter','none'); end
if ~isempty(ylab), ylabel(ylab,'interpreter','none'); end
if ndims(a) == 3 & ~isempty(clab)
  title({ clab ,titl },'interpreter','none');
else
  title(titl,'interpreter','none');
end

% ============================================================================


