function h=plot(a, varargin)
% h = plot(s, method) : plot iData object
%
%   @iData/plot function to plot data sets
%   This function plot the signal of the object as a function of the defined axes.
%   The plot is an errorbar plot for 1D data y=f(x), a surface mesh for 2D z=f(x,y),
%     and an isosurface for 3D c=f(x,y,z) data. Data specified as series of 
%     coordinate points are plotted using a plot3-type rendering. Further
%     dimensionalities are not handled.
%
%   The scatter3 rendering option is similar to plot3, but color points are set
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
%                 hide_errorbars is also valid not to plot error bars.
%               For 2D plots z=f(x,y), method is a string which may contain:
%                 surf, mesh, contour, contour3, surfc, surfl, contourf
%                 plot3, scatter3 (colored points), stem3, pcolor, waterfall
%               For 3D plots c=f(x,y,z), method is a string which may contain:
%                 plot3 (volume), scatter3 (colored points)
%                 waterfall, contour (set of coutour plots)
%                 surf, surf median, surf mean, surf half (isosurface)
%               The slice(a) method opens the interactive sliceomatic 3D viewer.
%
%               Global options for 2D and 3D plots: 
%                 flat, interp, faceted (for shading), view2, view3
%                 transparent, light, clabel, colorbar, shifted (overlayed 2D)
%               Global options for all plots: 
%                 axis tight, axis auto, hide_axes (compact layout)
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%                 whole (do not reduce object size for plotting)
%                 
% output: h: graphics object handles (cell/array)
% ex:     plot(iData(rand(10)), 'surfc interp transparent'); plot(iData(1:10), 'r-');
%         plot(iData(peaks));
%         [x,y,z,v]=flow; c=iData(x,y,z,v); plot(c,'surf');
%
% Contributed code (Matlab Central): 
%   fscatter3: Felix Morsdorf, Jan 2003, Remote Sensing Laboratory Zuerich
%   vol3d:     Joe Conti, 2004
%   sliceomatic: Eric Ludlam 2001-2008
%
% Version: $Revision: 1.75 $
% See also iData, interp1, interpn, ndgrid, plot, iData/setaxis, iData/getaxis
%          iData/xlabel, iData/ylabel, iData/zlabel, iData/clabel, iData/title
%          shading, lighting, surf, iData/slice

% private functions:
%   fscatter3: Felix Morsdorf, Jan 2003, Remote Sensing Laboratory Zuerich
%   vol3d:     Joe Conti, 2004

ih = ishold;
h  = {};
if nargin == 1, method=''; 
elseif length(varargin) == 1
  if ischar(varargin{1})
    method=varargin{1};
  elseif isa(varargin{1},'iData')
    a = [ a(:) ; varargin{1} ]; method='';
  end
else
  % split varvargin looking for chars
  method = '';
  index=1;
  while index <= length(varargin)  % parse input arguments and split with char/methods calls
    if ischar(varargin{index})
      method = varargin{index};
      h{end+1}=plot(a, method);
      a = []; method='';
      hold on
    elseif isa(varargin{index},'iData') 
      a = [a(:) ; varargin{index} ];
    end
    index=index+1;
  end
end

if isempty(a)
  if all(cellfun('length',h) == 1)
    h = cell2mat(h);
  end
  if ih == 1, hold on; else hold off; end
  return; 
end

% plot an array of objects
if length(a) > 1
  sum_max = 0;
  toremove='plot3 stem3 scatter3 mesh surf waterfall tight auto hide view2 view3 transparent';
  toremove=strread(toremove,'%s','delimiter',' ');
  this_method = method;
  for index=1:length(toremove)
    this_method = deblank(strrep(this_method, toremove{index},''));
  end
  for index=1:length(a(:))
    if ndims(a(index)) == 1 && isvector(a(index)) == 1 && ...
      isempty(getaxis(a(index),2)) && ...
      (~isempty(strfind(method, 'plot3'))      || ~isempty(strfind(method,'stem3')) ...
       || ~isempty(strfind(method,'scatter3')) || ~isempty(strfind(method, 'mesh')) ...
       || ~isempty(strfind(method,'surf') )    || ~isempty(strfind(method,'waterfall')))
      a(index) = setaxis(a(index), 2, index);
    end
    h{index} = plot(a(index), method);
    sum_max = sum_max+max(a(index))-min(a(index));
    if ndims(a(index)) == 1 && isempty(this_method)
      % change color of line
      colors = 'bgrcmk';
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
          set(h{index},'YData',z); 
        end
      end
    end
  end
  if all(cellfun('length',h) == 1)
    h = cell2mat(h);
  end
  h = reshape(h, size(a));
  if ih == 1, hold on; else hold off; end
  return
end

% plot a single object

% check if the object is not too large, else rebin accordingly
if prod(size(a)) > 1e6 && isempty(strfind(method,'whole'))
  iData_private_warning(mfilename, [ 'Object ' a.Tag ' is large (numel=' num2str(prod(size(a))) ...
    '.\n\tNow rebinning for display purposes with e.g. a=a(1:2:end, 1:2:end, ...).' ...
    '\n\tUse e.g plot(a, ''whole'') to plot the whole data set and be able to zoom tiny regions.' ]);
  a=iData_private_reduce(a);
end
zlab = '';
switch ndims(a) % handle different plotting methods depending on the iData dimensionality
case 0
  h=[]; 
  if ih == 1, hold on; else hold off; end
  return;
case 1  % vector type data (1 axis + signal) -> plot
  if size(a,1) ==1 && size(a,2) > 1
    a = transpose(a);
  end
  [x, xlab] = getaxis(a,1); x=double(x(:));
  [y, ylab] = getaxis(a,0); y=double(y(:));
  e         = get(a,'Error');   e=real(double(e)); e=e(:);
  m         = get(a,'Monitor'); m=real(double(m)); m=m(:);
  if not(all(m == 1 | m == 0)),
    e=genop(@rdivide,e,m); ylab = [ylab ' per monitor' ];
  end
  y=real(y);
  
  if isempty(method), method='b-'; end
  % handle side-by-side 1D plots
  if ~isempty(strfind(method,'plot3'))    | ~isempty(strfind(method,'stem3')) ...
   | ~isempty(strfind(method,'scatter3')) | ~isempty(strfind(method, 'mesh')) ...
   | ~isempty(strfind(method,'surf') )    | ~isempty(strfind(method,'waterfall'))
  	ax = getaxis(a,2);
  	if isempty(ax)
  		ax = 0;
    end
    if length(ax) == 1
    	ax = ax*ones(size(a));
    end
    % need to create this axis
    setalias(a, 'Axis_2', ax);
    setaxis(a, 2, 'Axis_2');
    h = plot(a, method);
    if ih == 1, hold on; else hold off; end
    return
  else 
    this_method=method;
    % clean up options that can be used together with linespec
    % tight, auto, tight, hide, view2, view3, transparent
    toremove='tight auto hide view2 view3 transparent';
    toremove=strread(toremove,'%s','delimiter',' ');
    for index=1:length(toremove)
      this_method = deblank(strrep(this_method, toremove{index},''));
    end
    if all(e == 0) || length(x) ~= length(e)
      if length(this_method)
        try
          h = plot(x,y, this_method);
        catch
          this_method=[];
        end
      end
      if ~length(this_method) h = plot(x,y); end
    else
      if length(this_method), 
        try
          h = errorbar(x,y,e,this_method);          
        catch
          this_method=[];
        end
      end
      if ~length(this_method) 
        h = errorbar(x,y,e); 
      end
      if ~isempty(strfind(this_method, 'hide_err')) || all(abs(e) >= abs(y) | e == 0)
        if length(h) == 1 && length(get(h,'Children') == 2)
          eh = get(h,'Children');
        else eh = h; 
        end
        set(eh(2), 'Visible','off');
      end
    end
  end
case 2  % surface type data (2 axes+signal) -> surf or plot3
  % check if a re-grid is needed
  if isvector(a) || (~isvector(a) && (~isempty(strfind(method,'plot3')) || ~isempty(strfind(method,'scatter3')) ))
    a = interp(a,'grid');
  end
  [x, xlab] = getaxis(a,2);
  [y, ylab] = getaxis(a,1);
  [z, zlab] = getaxis(a,0);
  m         = get(a,'Monitor');
  if not(all(m(:) == 1 | m(:) == 0)),
    zlab = [zlab ' per monitor' ];
  end
  x=real(double(x));
  y=real(double(y));
  z=real(double(z));
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
      if isempty(getaxis(a,3))
        [C,h]=contour(x,y,z);
      else
        c=getaxis(a,3); c=mean(c(:));
        Z(:,:,1)=z; Z(:,:,2)=z; Z(:,:,3)=z;
        h=contourslice(x,y,[c*0.999 c c*1.001],Z,[],[],c);
      end
    elseif (strfind(method,'surfc'))
      h    =surfc(x,y,z); % set(h,'Edgecolor','none');
    elseif (strfind(method,'surfl'))
      h    =surfl(x,y,z); set(h,'Edgecolor','none');
    elseif (strfind(method,'mesh'))
      h    =mesh(x,y,z);
    elseif ~isempty(strfind(method,'pcolor')) || ~isempty(strfind(method,'image'))
      h    =pcolor(x,y,z); set(h,'Edgecolor','none');
      if ~isempty(getaxis(a,3))
        c=getaxis(a,3); c=mean(c(:));
        zh= get(h,'ZData'); zh=ones(size(zh))*c;
        set(h,'ZData',zh);
      end
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
    elseif (strfind(method,'waterfall'))
      h=waterfall(x,y,z);
    else
      h=surf(x,y,z); set(h,'Edgecolor','none');
    end

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
    [x, xlab] = getaxis(a,2); x=double(x);
    [y, ylab] = getaxis(a,1); y=double(y);
    [z, zlab] = getaxis(a,3); z=double(z);
    [c, clab] = getaxis(a,0); c=double(c);
    m         = get(a,'Monitor');
    if not(all(m(:) == 1 | m(:) == 0)), clab = [clab ' per monitor' ]; end
    if isvector(a) == 3 || ~isempty(strfind(method, 'scatter3')) % plot3-like
      h=fscatter3(x(:),y(:),z(:),c(:));     % scatter3: require meshgrid
      view(3);
    else
      if ~isempty(strfind(method, 'plot3')) % vol3d: does not require meshgrid
        h = vol3d('cdata',c,'texture','3D','xdata',x,'ydata',y,'zdata',z);
        alphamap('vdown'); % make object transparent on borders and solid in center
        h = vol3d(h);
        h = h.handles;
      elseif ~isempty(strfind(method, 'waterfall')) || ~isempty(strfind(method, 'contour'))
        if ~isempty(strfind(method, ' y '))
          iy = linspace(min(y(:)), max(y(:)), 10);
          h = contourslice(x,y,z,c,[],iy,[]);
        elseif ~isempty(strfind(method, ' x '))
          ix = linspace(min(x(:)), max(x(:)), 10);
          h = contourslice(x,y,z,c,ix,[],[]);
        else
          iz = linspace(min(z(:)), max(z(:)), 10);
          h = contourslice(x,y,z,c,[],[],iz);
        end
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
        try
          if ~isempty(iso), 
            isosurface(x,y,z, c, iso);
          else 
            isosurface(x,y,z, c); 
          end
          h = findobj(gca,'type','patch');
        catch
          h = plot(a, 'scatter3');
          if ih == 1, hold on; else hold off; end
          return
        end
      end
    end
    zlabel(zlab);
  end
otherwise
  iData_private_warning(mfilename, [ 'plotting of ' num2str(ndims(a)) '-th dimensional data is not implemented from ' a.Tag '.\n\tUse sum or camproj to reduce dimensionality for plotting.' ]);
  h=[];
  if ih == 1, hold on; else hold off; end
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
if (strfind(method,'painters'))
	set(gcf,'Renderer','painters')
elseif (strfind(method,'zbuffer'))
	set(gcf,'Renderer','zbuffer')
end
if (strfind(method,'colorbar'))
  colorbar
end

% add a UIcontextMenu so that right-click gives info about the iData plot
T   = a.Title; if iscell(T), T=T{1}; end
T   = regexprep(T,'\s+',' '); % remove duplicated spaces
cmd = char(a.Command{end});
S   = a.Source;
[pS, fS, eS] = fileparts(S);
if length(pS) > 13, pS=[ '...' pS(end-10:end) ]; end
if length(fS) > 13, fS=[ '...' fS(end-10:end) ]; end
if ~isempty(pS), S = [ pS filesep ];
else             S = '';
end
S = [ S fS ];
if ~isempty(eS), S = [ S '.' eS ]; end
if length(cmd) > 23, cmd = [ cmd(1:20) '...' ]; end

% DisplayName and Label
d = '';
if ~isempty(a.Label) && ~isempty(a.DisplayName)
  g = cellstr(a.Label); g=deblank(g{1});
  if length(g) > 13, g = [ g(1:10) ]; end                 % Label/DisplayName
  d = [ d sprintf('%s', g) ];
  g = cellstr(a.DisplayName); g=deblank(g{1});
  if length(g) > 13, g = [ g(1:10) '...' ]; end           % 
  d = [ d sprintf('/%s', g) ];
elseif ~isempty(a.Label)
  g = cellstr(a.Label); g=deblank(g{1});
  if length(g) > 13, g = [ g(1:10) '...' ]; end           % Label
  d = [ d sprintf('%s', g) ];
elseif ~isempty(a.DisplayName)
  g = cellstr(a.DisplayName); g=deblank(g{1});
  if length(g) > 13, g = [ g(1:10) '...' ]; end           % DisplayName
  d = [ d sprintf('%s', g) ];
end

% --------------------- contextual menus ---------------------------------------
% create callback functions: handle callback for single error on/off (line
% or uimenu)
  function callback_toggle_error_gco(obj, event)
    if strcmp(get(obj,'type'),'uimenu'), obj=gco; end
    tmp_t = get(obj,'type');
    if strcmp(tmp_t,'hggroup')
      tmp_h=get(obj,'children'); % hggroup
    else
      return
    end
    % test if we have a line
    if length(tmp_h) ~= 2, return; end
    if ~strcmp(get(tmp_h(2),'type'),'line'), return; end
    % toggle visible
    if  strcmp(get(tmp_h(2),'visible'),'off'), tmp_v='on'; 
    else tmp_v='off'; end; 
    set(tmp_h(2),'visible',tmp_v);
  end

% create callback functions: handle callback for error on/off in axis frame
  function callback_toggle_error_gca(obj, event)
    % we scan all objects below the axis object, and toggle error bars
    hg = findobj(gca,'type','hggroup');
    for i=1:length(hg)
      callback_toggle_error_gco(hg(i))
    end
  end
  
% create callback functions: handle callback rotate 2d 3d
  function callback_rotate(obj, event)
    [tmp_a,tmp_e]=view; 
    if (tmp_a==0 & tmp_e==90) view(3); 
    else view(2); end; 
    lighting none;alpha(1);shading flat;rotate3d off;axis tight;
  end
  
% create callback functions: duplicate plot
  function callback_duplicate(obj, event)
    duplicate_cb.o =gco;
    duplicate_cb.g =gca;
    duplicate_cb.f =figure; 
    if strcmp(get(duplicate_cb.o,'type'),'axes')
      duplicate_cb.c=copyobj(duplicate_cb.g,gcf);
    else
      duplicate_cb.c=copyobj(duplicate_cb.o,gca);
    end
    set(gca,'position',[ 0.1 0.1 0.85 0.8]);
    set(gca,'XTickLabelMode','auto','XTickMode','auto');
    set(gca,'YTickLabelMode','auto','YTickMode','auto');
    set(gca,'ZTickLabelMode','auto','ZTickMode','auto');
    
    duplicate_cb.ud=get(duplicate_cb.g,'UserData'); 
    if ~isstruct(duplicate_cb.ud), return; end
    if iscellstr(duplicate_cb.ud.title) duplicate_cb.ud.title=duplicate_cb.ud.title{1}; end
    set(gcf,'Name', [ 'Copy of ' duplicate_cb.ud.title ]); 
    xlabel(duplicate_cb.ud.xlabel);ylabel(duplicate_cb.ud.ylabel); 
    title(duplicate_cb.ud.title);
    % create a simplified contextual menu
    uicm = uicontextmenu; 
    if ndims(a) >= 2 && ~isfield(ud,'contextual_2d')
      ud.contextual_2d = 1;
    end
    if isfield(duplicate_cb.ud,'contextual_2d') && duplicate_cb.ud.contextual_2d==1
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
    uimenu(uicm, 'Separator','on','Label','Toggle grid', 'Callback','grid');

    uimenu(uicm, 'Separator','on','Label', 'About iData', 'Callback',[ 'msgbox(''' version(iData)  sprintf('. Visit <http://ifit.mccode.org>') ''',''About iFit'',''help'')' ]);
    set(gca, 'UserData', ud);
set(gca, 'UIContextMenu', uicm);

   end
% create callback functions: about object
  function callback_about(obj, event)
    ud   = get(get(gco,'UIContextMenu'),'UserData');
    if ~isstruct(ud), return; end
    if iscellstr(ud.title) ud.title=ud.title{1}; end
    f = getframe(gcf);
    msgbox(  ud.properties, ...
             [ 'About: Figure ' num2str(gcf) ' ' ud.title ], ...
             'custom', f.cdata, get(gcf,'Colormap'));
  end
% create callback functions: about object
  function callback_commands(obj, event)
    ud   = get(get(gco,'UIContextMenu'),'UserData');
    if ~isstruct(ud), return; end
    titl = ud.title;
    if iscellstr(ud.title) ud.title=ud.title{1}; end
    listdlg('ListString', getfield(ud, 'commands'), ...
               'ListSize',[400 300],'Name', ud.title ,'PromptString', ud.name );
  end

% contextual menu for the single object being displayed
uicm = uicontextmenu; 
uimenu(uicm, 'Label', [ 'About ' a.Tag ': ' num2str(ndims(a)) 'D object ' mat2str(size(a)) ' ...' ], ...
  'Callback', @callback_about);
uimenu(uicm, 'Label',[ 'Duplicate "' T '" ...' ], 'Callback', @callback_duplicate);
if ndims(a) == 1
  uimenu(uicm, 'Label','Toggle error bars', 'Callback',@callback_toggle_error_gco);
end
uimenu(uicm, 'Separator','on', 'Label', [ 'Title: "' T '" ' d ], 'Callback', @callback_about);
if exist(a.Source,'file') && ~isdir(a.Source)
  uimenu(uicm, 'Label', [ 'Source: <' S '>' ], 'Callback',[ 'edit(''' a.Source ''')' ]);
else
  uimenu(uicm, 'Label', [ 'Source: <' S '>' ]);
end
uimenu(uicm, 'Label', [ 'Cmd: ' cmd ], 'Callback', @callback_commands);
uimenu(uicm, 'Label', [ 'User: ' a.User ], 'Callback', @callback_about);

% make up title string and Properties dialog content
properties={ [ 'Data ' a.Tag ': ' num2str(ndims(a)) 'D object ' mat2str(size(a)) ], ...
             [ 'Title: "' T '" ' d ], ...
             [ 'Source: ' a.Source ], ...
             [ 'Last command: ' cmd ]};

properties{end+1} = '[Rank]         [Value] [Description]';
uimenu(uicm, 'Separator','on', 'Label', '[Rank]         [Value] [Description]');
for index=0:length(getaxis(a))
  [v, l] = getaxis(a, num2str(index));
  if length(l) > 20, l = [l(1:18) '...' ]; end 
  x      = getaxis(a, index);
  m      = get(a, 'Monitor');
  if index==0 & not(all(m==1 | m==0))
    t = sprintf('%6i %15s  %s [%g:%g] (per monitor=%g)', index, v, l, full(min(x(:))), full(max(x(:))), mean(m(:)));
  else
    try
      [s, f] = std(a, index);
      t = sprintf('%6i %15s  %s [%g:%g] <%g +/- %g>', index, v, l, full(min(x(:))), full(max(x(:))), f, s);
    catch
      t = sprintf('%6i %15s  %s [%g:%g]', index, v, l, full(min(x(:))), full(max(x(:))));
    end
  end
  properties{end+1} = t;
  uimenu(uicm, 'Label', t);
end
titl ={ T ; [ a.Tag ' <' S '>' ]};
if length(T) > 23, T=[ T(1:20) '...' ]; end
if length(S)+length(d) < 30,
  d = [ d ' ' T ];
end
try
  if ~isempty(d)
    set(h, 'DisplayName', [ d ' ' a.Tag ' <' S '>']);
  else
    set(h, 'DisplayName', [ T a.Tag ' <' S '>' ]);
  end
catch
end
uimenu(uicm, 'Separator','on','Label', 'About iData', 'Callback',[ 'msgbox(''' version(iData)  sprintf('. Visit <http://ifit.mccode.org>') ''',''About iFit'',''help'')' ]);
% attach contexual menu to plot with UserData storage
ud.properties=properties;
ud.xlabel = xlab;
ud.ylabel = ylab;
ud.zlabel = zlab;
ud.title  = titl;
ud.name   = char(a);
ud.commands = commandhistory(a);
set(uicm,'UserData', ud);
set(h,   'UIContextMenu', uicm); 

% contextual menu for the axis frame
if ~isempty(get(gca, 'UserData'))
  ud = get(gca, 'UserData');
end
uicm = uicontextmenu;
uimenu(uicm, 'Label',[ 'Duplicate this view...' ], 'Callback', @callback_duplicate);;
if ndims(a) == 1 && ~isfield(ud,'contextual_1d')
  ud.contextual_1d = 1;
end
if isfield(ud,'contextual_1d') && ud.contextual_1d==1
  uimenu(uicm, 'Label','Toggle All Error Bars', 'Callback', @callback_toggle_error_gca);
end
uimenu(uicm, 'Label','Toggle grid', 'Callback','grid');
if ndims(a) >= 2 && ~isfield(ud,'contextual_2d')
  ud.contextual_2d = 1;
end
if isfield(ud,'contextual_2d') && ud.contextual_2d==1
  uimenu(uicm, 'Label','Reset Flat/3D View', 'Callback',@callback_rotate);
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
  set(gcf, 'KeyPressFcn', @(src,evnt) eval('if lower(evnt.Character)==''r'', lighting none;alpha(1);shading flat;axis tight;rotate3d off; end') );
  if ndims(a) >= 2
    uimenu(uicm, 'Label', 'Rotate', 'Callback','rotate3d on');
  end
end
uimenu(uicm, 'Separator','on','Label', 'About iData', 'Callback',[ 'msgbox(''' version(iData) sprintf('. Visit <http://ifit.mccode.org>') ''',''About iFit'',''help'')' ]);

set(gca, 'UserData', ud);
set(gca, 'UIContextMenu', uicm);
try
  set(h,   'Tag',  [ mfilename '_' a.Tag ]);
end
set(gcf, 'Name', char(a));

% labels
if (strfind(method,'hide_ax'))
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
  if ~isempty(d)
    titl{1} = [ titl{1} ' ''' d '''' ]; 
  end
  title(titl,'interpreter','none');
end

if ih == 1, hold on; else hold off; end

end
% ============================================================================


