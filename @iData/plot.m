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
%   The scatter3 rendering option is similar to plot3, but colors points 
%   according to the signal intensity. The 'plot3' option for 3D (volume) objects
%   uses a semi-transparent volume rendering, whereas the default plot uses
%   an iso-surface on the median signal.
%
%   As mentioned in the iData axis definition (see iData/setaxis), the 'X' axis
%     refers to the 2nd dimension (along columns), whereas the 'Y' axis refers to
%     the first dimension (along rows).
%
% input:  s: object or array (iData)
%         method: optional type of plot to render
%
%               For 1D plots y=f(x), method is a string to specify color/symbol.
%               For 2D plots z=f(x,y), method is a string which may contain:
%                 surf, mesh, contour, contour3, surfc, surfl, contourf
%                 plot3, scatter3 (colored points), stem3, pcolor
%               For 3D plots c=f(x,y,z), method is a string which may contain:
%                 plot3 (volume), scatter3 (colored points)
%                 surf, surf median, surf mean, surf half
%               The slice(a) method opens the interactive sliceomatic 3D viewer.
%
%               Global options for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel, colorbar, shifted (overlayed 2D)
%                 axis tight, axis auto, view2, view3, hide_axes (compact layout)
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%                 
% output: h: graphics object handles (cell)
% ex:     plot(iData(rand(10)), 'surfc interp transparent'); plot(iData(1:10), 'r-');
%         plot(iData(peaks));
%         [x,y,z,v]=flow; c=iData(x,y,z,v); plot(c,'surf');
%
% Contributed code (Matlab Central): 
%   fscatter3: Felix Morsdorf, Jan 2003, Remote Sensing Laboratory Zuerich
%   vol3d:     Joe Conti, 2004
%   sliceomatic: Eric Ludlam 2001-2008
%
% Version: $Revision: 1.49 $
% See also iData, interp1, interpn, ndgrid, plot, iData/setaxis, iData/getaxis
%          iData/xlabel, iData/ylabel, iData/zlabel, iData/clabel, iData/title
%          shading, lighting, surf, iData/slice

% private functions:
%   fscatter3: Felix Morsdorf, Jan 2003, Remote Sensing Laboratory Zuerich
%   vol3d:     Joe Conti, 2004
if nargin == 1, method=''; end
if length(a) > 1
  h = cell(size(a));
  ih = ishold;
  sum_max = 0;
  for index=1:length(a(:))
    h{index} = plot(a(index), method);
    sum_max = sum_max+max(a(index))-min(a(index));
    if ndims(a(index)) == 1 && isempty(method)
      % change color of line
      colors = 'bgrcm';
      set(h{index}, 'color', colors(1+mod(index, length(colors))));
    end
    hold on
  end
  % re-arrange if this is a 2D overlay
  for index=1:length(a(:))
    if length(h{index}) == 1 && ~isempty(strfind(method, 'shifted'))
      if ndims(a(index)) ~= 1
        try
          z= get(h{index},'ZData');
          if all(z == 0), use_cdata=1; z= get(h{index},'CData');
          else use_cdata=0; end
          z = z-min(z(:));
          z = z+sum_max*index/length(a(:));
          if use_cdata==0, set(h{index},'ZData',z);
          else set(h{index},'CData',z); end
        end
      else
        try
          z= get(h{index},'YData');
          z = z-min(z(:));
          z = z+sum_max*index/length(a(:));
          set(h{index},'YData',z); end
        end
      end
    end
  end
  h = reshape(h, size(a));
  if ih == 1, hold on; else hold off; end
  return
end

% check if the object is not too large, else rebin accordingly
if prod(size(a)) > 1e6
  iData_private_warning(mfilename, [ 'Object ' a.Tag ' is too large (numel=' num2str(prod(size(a))) ...
    '.\n\tYou should rebin with e.g. a=a(1:2:end, 1:2:end, ...).' ]);
end

switch ndims(a) % handle different plotting methods depending on the iData dimensionality
case 0
  h=[]; return;
case 1  % vector type data (1 axis + signal) -> plot
  if size(a,1) ==1 && size(a,2) > 1
    a = transpose(a);
  end
  [x, xlab] = getaxis(a,1); x=x(:);
  [y, ylab] = getaxis(a,0); y=y(:);
  e         = get(a,'Error');   e=real(e); e=e(:);
  m         = get(a,'Monitor'); m=real(m); m=m(:);
  if not(all(m == 1 | m == 0)),
    e=genop(@rdivide,e,m); ylab = [ylab ' per monitor' ];
  end
  y=real(y);
  
  if isempty(method), method='b-'; end
  % handle side-by-side 1D plots
  if ~isempty(strfind(method, 'plot3')) | ~isempty(strfind(method,'stem3')) ...
   | ~isempty(strfind(method,'scatter3')) | ~isempty(strfind(method, 'mesh')) | ~isempty(strfind(method,'surf') )
  	if isempty(getaxis(a,2))
  		ax = 0;
    else
    	ax = getaxis(a, 2);
    end
    if length(ax) == 1
    	ax = ax*ones(size(a));
    end
    % need to create this axis
    setalias(a, 'Axis_2', ax);
    setaxis(a, 2, 'Axis_2');
    h = plot(a, method);
    return
  else 
    this_method=method;
    if all(e == 0) | length(x) ~= length(e)
      if length(this_method)
        try
          h = plot(x,y, method);
        catch
          this_method=[];
        end
      end
      if ~length(this_method) h = plot(x,y); end
    else
      if length(method), 
        try
          h = errorbar(x,y,e,method);
        catch
          this_method=[];
        end
      end
      if ~length(this_method) h = errorbar(x,y,e); end
    end
  end
  set(h, 'Tag', a.Tag);
case 2  % surface type data (2 axes+signal) -> surf or plot3
  % check if a re-grid is needed
  if isvector(a) || (~isvector(a) && (~isempty(strfind(method,'plot3')) || ~isempty(strfind(method,'scatter3')) ))
    a = interp(a,'grid');
  end
  [x, xlab] = getaxis(a,2);
  [y, ylab] = getaxis(a,1);
  [z, zlab] = getaxis(a,0);
  m         = get(a,'Monitor');
  if not(all(m == 1 | m == 0)),
    zlab = [zlab ' per monitor' ];
  end
  x=real(x);
  y=real(y);
  z=real(z);
  if isvector(a) % plot3/fscatter3
    if (strfind(method,'plot3'))
    	method = strrep(method,'plot3','');
    	method = strrep(method,' ','');
      if length(method), h = plot3(x,y,z, method);
      else h = plot3(x,y,z); end
    else
      h=fscatter3(x,y,z,z); view(3);
    end
  else                % surf and similar stuff
    C = [];
    if isvector(x) & isvector(y),
      z = z;
    end
    if (strfind(method,'contour3'))
      [C,h]=contour3(x,y,z);
    elseif (strfind(method,'contourf'))
      [C,h]=contourf(x,y,z);
    elseif (strfind(method,'contour'))
      [C,h]=contour(x,y,z);
    elseif (strfind(method,'surfc'))
      h    =surfc(x,y,z); % set(h,'Edgecolor','none');
    elseif (strfind(method,'surfl'))
      h    =surfl(x,y,z); set(h,'Edgecolor','none');
    elseif (strfind(method,'mesh'))
      h    =mesh(x,y,z);
    elseif ~isempty(strfind(method,'pcolor')) || ~isempty(strfind(method,'image'))
      h    =pcolor(x,y,z); set(h,'Edgecolor','none');
    elseif (strfind(method,'stem3'))
    	method = strrep(method,'stem3','');
    	method = strrep(method,' ','');
      if length(method), h = stem3(x,y,z, method);
      else h = stem3(x,y,z); end
    elseif (strfind(method,'plot3'))
      a = interp(a,'grid');
    	method = strrep(method,'plot3','');
    	method = strrep(method,' ','');
      if length(method), h = plot3(x(:),y(:),z(:), method);
      else h = plot3(x,y,z); end
    elseif (strfind(method,'scatter3'))
      h=fscatter3(x(:),y(:),z(:),z(:));
    else
      h=surf(x,y,z); set(h,'Edgecolor','none');
    end
    set(h, 'Tag', a.Tag);

    if ~isempty(C) & strfind(method,'clabel')
      clabel(C,h);
    end
  end
  zlabel(zlab);
case 3  % 3d data sets: volumes
  % first test if this is an image
  if isfield(a.Data,'cdata')
    h=image(a.Data.cdata);
    xlab=''; ylab=''; clab='';
  else
    % check if a rebining on a grid is required
    if ~isvector(a) && isempty(strfind(method, 'plot3'))
      a = interp(a,'grid'); % make sure we get a grid
    end
    [x, xlab] = getaxis(a,2);
    [y, ylab] = getaxis(a,1);
    [z, zlab] = getaxis(a,3);
    [c, clab] = getaxis(a,0);
    m         = get(a,'Monitor');
    if not(all(m == 1 | m == 0)), clab = [clab ' per monitor' ]; end
    if isvector(a) == 3 || ~isempty(strfind(method, 'scatter3')) % plot3-like
      h=fscatter3(x(:),y(:),z(:),c(:));     % scatter3: require meshgrid
      view(3);
    else
      if ~isempty(strfind(method, 'plot3')) % vol3d: does not require meshgrid
        h = vol3d('cdata',c,'texture','3D','xdata',x,'ydata',y,'zdata',z);
        alphamap('vdown'); % make object transparent on borders and solid in center
        h = vol3d(h);
        h = h.handles
      elseif ~isempty(strfind(method, 'slice')) % sliceomatic
        slice(a); h=[];
      else
        % isosurface: require meshgrid
        if ~isempty(strfind(method, 'mean'))
          iso = mean(c(:));
        elseif ~isempty(strfind(method, 'half'))
          iso = (min(c(:))+max(c(:)))/2;
        elseif ~isempty(strfind(method, 'median'))
          iso = median(c(:));
        else
          iso = [];
        end
        if ~isempty(iso), 
          isosurface(x,y,z, c, iso);
        else 
          isosurface(x,y,z, c); 
        end
        h = findobj(gca,'type','patch');
      end
    end
    try
      set(h, 'Tag', a.Tag);
    end
    zlabel(zlab);
  end
otherwise
  iData_private_warning(mfilename, [ 'plotting of ' num2str(ndims(a)) '-th dimensional data is not implemented from ' a.Tag '.\n\tUse sum or camproj to reduce dimensionality for plotting.' ]);
  h=[];
  return
end % switch

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
if (strfind(method,'view2'))
  view(2);
end
if (strfind(method,'view3'))
  view(3);
end
if (strfind(method,'tight'))
  axis tight
end
if (strfind(method,'auto'))
  axis auto
end
if findstr(method,'painters')
	set(gcf,'Renderer','painters')
elseif findstr(method,'zbuffer')
	set(gcf,'Renderer','zbuffer')
end
if (strfind(method,'colorbar'))
  colorbar
end

% add a UIcontextMenu so that right-click gives info about the iData plot
% also make up title string and Properties dialog content
T=a.Title; if iscell(T), T=T{1}; end
T=regexprep(T,'\s+',' '); % remove duplicated spaces
cmd = char(a.Command{end});
properties={ [ 'Data ' a.Tag ': ' num2str(ndims(a)) 'D object ' mat2str(size(a)) ], ...
             [ 'Title: "' T '"' ], ...
             [ 'Source: ' a.Source ], ...
             [ 'Last command: ' cmd ]};
if ~isempty(a.Label)
  properties{end+1} = [ 'Label: ' a.Label ];
end
if length(getaxis(a))
  properties{end+1} = '[Rank]         [Value] [Description]';
  for index=0:length(getaxis(a))
    [v, l] = getaxis(a, num2str(index));
    x      = getaxis(a, index);
    m      = get(a, 'Monitor');
    if index==0 & not(all(m==1 | m==0))
      properties{end+1} = sprintf('%6i %15s  %s [%g:%g] (per monitor)', index, v, l, min(x(:)), max(x(:)));
    else
      properties{end+1} = sprintf('%6i %15s  %s [%g:%g]', index, v, l, min(x(:)), max(x(:)));
    end
  end
end
S=a.Source;
[dummy, fS] = fileparts(S);
titl ={ T ; [ a.Tag ' <' fS '>' ]};
if length(T) > 23, T=[ T(1:20) '...' ]; end
try
  if ~isempty(a.DisplayName)
    set(h, 'DisplayName', [ a.DisplayName ' ' T ]);
  else
    set(h, 'DisplayName', [ T a.Tag ' <' fS '>' ]);
  end
catch
end
if length(S) > 23, S=[ '...' S(end-20:end) ]; end
if length(cmd) > 23, cmd = [ cmd(1:20) '...' ]; end

uicm = uicontextmenu; 
set(uicm,'UserData', properties); 
uimenu(uicm, 'Label', [ 'About ' a.Tag ': ' num2str(ndims(a)) 'D object ' mat2str(size(a)) ' ...' ], ...
  'Callback', [ 'msgbox(get(get(gco,''UIContextMenu''),''UserData''), ''About: Figure ' num2str(gcf) ' ' T ' <' S '>'',''help'');' ] );
uimenu(uicm, 'Label',[ 'Duplicate ' T ' ...' ], 'Callback',[ 'g=gca; f=figure; c=copyobj(g,f); set(c,''position'',[ 0.1 0.1 0.9 0.8]); set(f,''Name'',''Copy of ' char(a) ''');' ]);
uimenu(uicm, 'Separator','on', 'Label', [ 'Title: "' T '"' ]);
uimenu(uicm, 'Label', [ 'Source: <' S '>' ]);
uimenu(uicm, 'Label', [ 'Cmd: ' cmd ]);
uimenu(uicm, 'Label', [ 'User: ' a.User ]);
uimenu(uicm, 'Separator','on','Label','Toggle grid', 'Callback','grid');
if ndims(a) == 1
  uimenu(uicm, 'Label','Toggle error bars', 'Callback','tmp_h=get(gco,''children''); if strcmp(get(tmp_h(2),''visible''),''off''), tmp_v=''on''; else tmp_v=''off''; end; set(tmp_h(2),''visible'',tmp_v); clear tmp_h tmp_v');
end
if ndims(a) >= 2
  uimenu(uicm, 'Label','Reset Flat/3D View', 'Callback','[tmp_a,tmp_e]=view; if (tmp_a==0 & tmp_e==90) view(3); else view(2); end; clear tmp_a, tmp_e; lighting none;alpha(1);shading flat;rotate3d off;axis tight;');
  uimenu(uicm, 'Label','Smooth View','Callback', 'shading interp;');
  uimenu(uicm, 'Label','Add Light','Callback', 'light;lighting phong;');
  uimenu(uicm, 'Label','Transparency','Callback', 'alpha(0.7);');
  uimenu(uicm, 'Label','Linear/Log scale','Callback', 'if strcmp(get(gca,''zscale''),''linear'')  set(gca,''zscale'',''log''); else set(gca,''zscale'',''linear''); end');
  uimenu(uicm, 'Label','Toggle Perspective','Callback', 'if strcmp(get(gca,''Projection''),''orthographic'')  set(gca,''Projection'',''perspective''); else set(gca,''Projection'',''orthographic''); end');
else
  uimenu(uicm, 'Label','Reset View', 'Callback','view(2);lighting none;alpha(1);shading flat;axis tight;rotate3d off;');
  uimenu(uicm, 'Label','Linear/Log scale','Callback', 'if strcmp(get(gca,''yscale''),''linear'')  set(gca,''yscale'',''log''); else set(gca,''yscale'',''linear''); end');
end
% add rotate/pan/zoom tools in case java machine is not started
if ~usejava('jvm')
  uimenu(uicm, 'Separator','on','Label','Zoom on/off', 'Callback','zoom');
  uimenu(uicm, 'Label','Pan (move)', 'Callback','pan');
  set(gcf, 'KeyPressFcn', @(src,evnt) eval('if lower(evnt.Character)==''r'', lighting none;alpha(1);shading flat;axis tight;rotate3d off;; end') );
  if ndims(a) >= 2
    uimenu(uicm, 'Label', 'Rotate', 'Callback','rotate3d on');
  end
end
uimenu(uicm, 'Label', 'About iData', 'Callback',[ 'msgbox(''' version(iData) ''')' ]);
% attach contexual menu to plot
set(h,   'UIContextMenu', uicm); 
set(gca, 'UIContextMenu', uicm);
set(h,   'Tag',  char(a));
set(gcf, 'Name', char(a));

% labels
if findstr(method,'hide_axes')
  % set(gca,'visible','off'); 
  set(gca,'XTickLabel',[],'XTick',[]); set(gca,'YTickLabel',[],'YTick',[]); set(gca,'ZTickLabel',[],'ZTick',[])
  xlabel(' '); ylabel(' '); zlabel(' ');
  title(a.Title,'interpreter','none');
else
  if ~isempty(xlab), xlabel(xlab,'interpreter','none'); end
  if ~isempty(ylab), ylabel(ylab,'interpreter','none'); end
  if ndims(a) == 3 & ~isempty(clab)
    titl = { clab , titl{:} };
  end
  title(titl,'interpreter','none');
end

% ============================================================================


