function h=plot(a, varargin)
% h = plot(s, method, ...) : plot iData object
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
%  Type <a href="matlab:doc(iData,'Plot')">doc(iData,'Plot')</a> to access the iFit/Plot Documentation.
%
% input:  s: object or array (iData)
%         method: optional type of plot to render
%
%               For 1D plots y=f(x), method is a string to specify color/symbol.
%                 hide_errorbars is also valid not to plot error bars.
%                 To plot a set of 1D objects side by side, specify a 2D plot 
%                 option such as 'surf' or 'plot3'.
%               For 2D plots z=f(x,y), method is a string which may contain:
%                 surf, mesh, contour, contour3, surfc, surfl, contourf
%                 plot3, scatter3 (colored points), stem3, pcolor, waterfall
%               For 3D plots c=f(x,y,z), method is a string which may contain:
%                 plot3 (volume), scatter3 (colored points, supports 'scatter3 filled' and 'scatter3 bubble')
%                 waterfall (supports 'waterfall x' 'y' and 'z'), contour (set of coutour plots)
%                 surf, surf median, surf mean, surf half (isosurface)
%               The slice(a) method opens the interactive sliceomatic 3D viewer.
%
%               Global options for 2D and 3D plots: 
%                 flat, interp, faceted (for shading), view2, view3
%                 transparent, light, clabel, colorbar, shifted (overlayed 2D)
%               Global options for all plots: 
%                 axis tight, axis auto, hide_axes (compact layout)
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%                 opengl (faster for large data sets)
%                 whole or full (do not reduce large object size for plotting)
%                 figure (open 
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
% Version: $Revision: 1.97 $
% See also iData, interp1, interpn, ndgrid, plot, iData/setaxis, iData/getaxis
%          iData/xlabel, iData/ylabel, iData/zlabel, iData/clabel, iData/title
%          shading, lighting, surf, iData/slice

% private functions:
%   fscatter3: Felix Morsdorf, Jan 2003, Remote Sensing Laboratory Zuerich
%   vol3d:     Joe Conti, 2004

ih = ishold;
h  = [];

% analyze input arguments
if nargin == 1, method=''; 
elseif length(varargin) == 1
  if ischar(varargin{1})
    method=varargin{1};
  elseif isa(varargin{1},'iData')
    b = varargin{1};
    a = [ a(:) ; b(:) ]; method='';
  end
else
  % split varargin looking for chars
  method = '';
  index=1;
  while index <= length(varargin)  % parse input arguments and split with char/methods calls
    if ischar(varargin{index})
      method = varargin{index};
      h =[ h plot(a, method) ];
      a = []; method='';
      hold on
    elseif isa(varargin{index},'iData') 
      b = varargin{index};
      a = [ a(:) ; b(:) ];
    end
    index=index+1;
  end
end
clear varargin

if isempty(a)
  if ih == 1, hold on; else hold off; end
  return; 
end

% clean method string from the plot type and supported options not to be passed to matlab plot commands
if ischar(method)
  toremove='plot3 stem3 scatter3 mesh surf waterfall tight auto hide view2 view3 transparent axis hide_err contour contour3 surfc surfl contourf pcolor median mean half slice flat interp faceted light clabel colorbar shifted hide_axes painters zbuffer whole full';
  toremove=strread(toremove,'%s','delimiter',' ');
  this_method = method;
  for index=1:length(toremove)
    [d1,d2,d3,d4]= regexp(this_method,[ '\<' toremove{index} ]);
    if isempty(d1), continue; end
    next_space   = find(this_method(d1:end) == ' ');
    if isempty(next_space), d2=length(this_method);
    elseif length(next_space) >= 1, d2=d1+next_space(1)-1; end
    this_method(d1:d2)='';
  end
else
  this_method = method;
end

% plot an array of objects
if numel(a) > 1
  iData_private_warning('enter', mfilename);
  sum_max = 0;
  % plot objects in the same axis frame
  % set error bar uniformly along objects
  common_error_bar='undefined'; % will set the value to 0/1 when 1D found
  for index=1:numel(a(:))
    if isempty(a(index)), h{index} = []; continue; end
    if ndims(a(index)) == 1 && isvector(a(index)) == 1 && ...
      isempty(getaxis(a(index),2)) && ...
      (~isempty(strfind(method, 'plot3'))      || ~isempty(strfind(method, 'stem3')) ...
       || ~isempty(strfind(method,'scatter3')) || ~isempty(strfind(method, 'mesh')) ...
       || ~isempty(strfind(method,'surf') )    || ~isempty(strfind(method, 'waterfall')))
      a(index) = setaxis(a(index), 2, index);
    end
    h{index} = plot(a(index), method);
    if ndims(a(index)) == 1
      if length(h{index}) == 1 && length(get(h{index},'Children') == 2)
        eh = get(h{index},'Children');
      else eh = h{index}; 
      end
      if length(eh) > 1
        if strcmp(common_error_bar, 'undefined')
          common_error_bar = get(eh(2), 'Visible');
        else
          set(eh(2), 'Visible',common_error_bar);
        end
      end
    end
    s = getaxis(a(index), 0);
    sum_max = sum_max+max(s(:))-min(s(:));
    this_h = get(h{index},'Type'); if iscell(this_h), this_h=this_h{1}; end
    if ndims(a(index)) == 1 && isempty(this_method) ...
    && any(strcmp(this_h,{'line','hggroup'})) && isempty(strfind(method,'scatter3'))
      % change color of line
      colors = 'bgrcmk';
      set(h{index}, 'color', colors(1+mod(index, length(colors))));

    end
    hold on
  end % for
  
  % re-arrange if this is a 2D overlay (shifted)
  if all(cellfun('length',h) <= 1)
    h = cell2mat(h);
  end
  for index=1:numel(h)
    if length(h(index)) == 1 && ~isempty(strfind(method, 'shifted'))
      if ndims(a(index)) ~= 1
        try
          z= get(h(index),'ZData'); 
          c= get(h(index),'CData');
          if all(z(:) == 0)
               use_cdata=1; z= c;
          else use_cdata=0; 
          clear c
          end
          z = z-min(z(:));
          z = z+sum_max*index/numel(a);
          if use_cdata==0, 
               set(h(index),'ZData',z);
          else set(h(index),'CData',z); 
          end
        end
      else
        try
          z= get(h(index),'YData');
          z = z-min(z(:));
          z = z+sum_max*index/length(a(:));
          set(h(index),'YData',z); 
        end
      end
    end
  end
  
  if ih == 1, hold on; else hold off; end
  
  iData_private_warning('exit', mfilename);
  return
end % plot array

% plot a single object
method = lower(method);

% check if the object is not too large, else rebin accordingly
if prod(size(a)) > 1e6 
  if isempty([ strfind(method,'whole') strfind(method,'full') ])
    iData_private_warning(mfilename, [ 'Object ' a.Tag ' is large (numel=' num2str(prod(size(a))) ...
      '.\n\tNow rebinning for display purposes with e.g. a=a(1:2:end, 1:2:end, ...).' ...
      '\n\tUse e.g plot(a, ''whole'') to plot the whole data set and be able to zoom tiny regions.' ]);
    a=iData_private_reduce(a);
  else
    method = [ method ' opengl' ];
  end
end
zlab = '';

% possibly select Rendered prior to start plotting
if (strfind(method,'opengl'))   % faster for large data sets
	set(gcf,'Renderer','OpenGL')
elseif (strfind(method,'painters'))
	set(gcf,'Renderer','painters')
elseif (strfind(method,'zbuffer'))
	set(gcf,'Renderer','zbuffer');
end

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
   | ~isempty(strfind(method,'scatter3')) | ~isempty(strfind(method,'mesh')) ...
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
      if ~isempty(strfind(method, 'hide_err')) || all(abs(e) >= abs(y) | e == 0)
        if length(h) == 1 && length(get(h,'Children') == 2)
          eh = get(h,'Children');
        else eh = h; 
        end
        if length(eh) > 1, set(eh(2), 'Visible','off'); end
      end
    end
  end
  clear x y e m
case 2  % surface type data (2 axes+signal) -> surf or plot3
  % check if a re-grid is needed
  if isvector(a)
    a_is_vector = 1; % plot as lines
  elseif (~isvector(a) && (~isempty(strfind(method,'plot3')) || ~isempty(strfind(method,'scatter3')) ))
    a = interp(a,'grid');
    a_is_vector = 1; % plot as lines, even after re-sampling (requested explicitly)
  else
    a_is_vector = 0;
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
  if a_is_vector % plot3/fscatter3
    if (strfind(method,'scatter3'))
      h=fscatter3(x(:),y(:),z(:),z(:),this_method); view(3);
    else
      if length(method), h = plot3(x,y,z, this_method);
      else h = plot3(x,y,z); end
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
      if length(method), h = stem3(x,y,z, this_method);
      else h = stem3(x,y,z); end
    elseif (strfind(method,'plot3'))
      a = interp(a,'grid');
      if length(method), h = plot3(x(:),y(:),z(:), this_method);
      else h = plot3(x,y,z); end
    elseif (strfind(method,'scatter3'))
      h=fscatter3(x(:),y(:),z(:),z(:),this_method);
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
  clear x y z m
case 3  % 3d data sets: volumes
  % first test if this is an image
  if isfield(a.Data,'cdata')
    h=image(a.Data.cdata);
    xlab=''; ylab=''; clab='';
  else
    % check if a rebining on a grid is required
    if ~isvector(a) && isempty(strfind(method, 'plot3')) && isempty(strfind(method, 'scatter3')) 
    whos a
      a = interp(a,1); % make sure we get a grid
    end
    [x, xlab] = getaxis(a,2); x=double(x);
    [y, ylab] = getaxis(a,1); y=double(y);
    [z, zlab] = getaxis(a,3); z=double(z);
    [c, clab] = getaxis(a,0); c=double(c);
    m         = get(a,'Monitor');
    if not(all(m(:) == 1 | m(:) == 0)), clab = [clab ' per monitor' ]; end
    if isvector(a) >= 3 || ~isempty(strfind(method, 'scatter3')) % plot3-like
      if ~isempty(strfind(method, 'scatter3'))
        h=fscatter3(x(:),y(:),z(:),c(:), this_method);     % scatter3: may require meshgrid
      else
        h=plot3(x(:),y(:),z(:), this_method);
      end
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
    clear x y z c m
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
if (strfind(method,'colorbar'))
  colorbar
end

% add a UIcontextMenu so that right-click gives info about the iData plot
T   = a.Title; if ~ischar(T), T=char(T); end
if ~isvector(T), T=transpose(T); T=T(:)'; end
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
  if length(g) > 23, g = [ g(1:20) '...' ]; end           % Label
  d = [ d sprintf('%s', g) ];
elseif ~isempty(a.DisplayName)
  g = cellstr(a.DisplayName); g=deblank(g{1});
  if length(g) > 23, g = [ g(1:20) '...' ]; end           % DisplayName
  d = [ d sprintf('%s', g) ];
end

T_char = char(T);
titl ={ T ; [ a.Tag ' <' S '>' ]};
if length(T) > 23, T=[ T(1:20) '...' ]; end
if length(S)+length(d) < 30,
  d = [ d ' ' T ];
end
try
  if ~isempty(d)
    set(h, 'DisplayName', [ d ]);
  else
    set(h, 'DisplayName', [ T a.Tag ' <' S '>' ]);
  end
catch
end

% contextual menu for the single object being displayed
% internal functions must be avoided as it uses LOTS of memory
uicm = uicontextmenu; 
% menu About
uimenu(uicm, 'Label', [ 'About ' a.Tag ': ' num2str(ndims(a)) 'D object ' mat2str(size(a)) ' ...' ], ...
  'Callback', [ 'msgbox(getfield(get(get(gco,''UIContextMenu''),''UserData''),''properties''),' ...
                '''About: Figure ' num2str(gcf) ' ' T ' <' S '>'',' ...
                '''custom'',getfield(getframe(gcf),''cdata''), get(gcf,''Colormap''));' ] );

% menu Toggle error bars (1D object)
if ndims(a) == 1
  uimenu(uicm, 'Label','Toggle Error Bars', 'Callback', [...
   'tmp_h=get(gco,''children'');'...
   'if strcmp(get(tmp_h(2),''visible''),''off''), tmp_v=''on''; else tmp_v=''off''; end;' ...
   'set(tmp_h(2),''visible'',tmp_v); clear tmp_h tmp_v' ]);
end
uimenu(uicm, 'Separator','on', 'Label', [ 'Title: "' T '" ' d ]);
if exist(a.Source,'file') && ~isdir(a.Source)
  uimenu(uicm, 'Label', [ 'Source: <' S '>' ], 'Callback',[ 'edit(''' a.Source ''')' ]);
else
  uimenu(uicm, 'Label', [ 'Source: <' S '>' ]);
end
if ~isempty(d)
  uimenu(uicm, 'Label', [ 'Label: ' d ]);
end
% menu List of Commands (history)
uimenu(uicm, 'Label', [ 'Cmd: ' cmd ' ...' ], 'Callback', [ ...
 'tmp_ud = get(get(gco,''UIContextMenu''),''UserData'');' ...
 'listdlg(''ListString'', getfield(tmp_ud, ''commands''),' ...
  '''ListSize'',[400 300],''Name'', tmp_ud.title ,''PromptString'', tmp_ud.name );clear tmp_ud' ])
uimenu(uicm, 'Label', [ 'User: ' a.User ]);

% make up title string and Properties dialog content
properties={ [ 'Data ' a.Tag ': ' num2str(ndims(a)) 'D object ' mat2str(size(a)) ], ...
             [ 'Title: "' T_char '" ' d ], ...
             [ 'Source: ' a.Source ], ...
             [ 'Last command: ' cmd ]};

properties{end+1} = '[Rank]         [Value] [Description]';
uimenu(uicm, 'Separator','on', 'Label', '[Rank]         [Value] [Description]');
for index=0:length(getaxis(a))
  [v, l] = getaxis(a, num2str(index));
  if length(l) > 20, l = [l(1:18) '...' ]; end 
  x      = getaxis(a, index);
  m      = get(a, 'Monitor');
  if length(x) == 1
    minmaxstd = sprintf('[%g]', full(x));
  elseif isvector(x)
    minmaxstd = sprintf('[%g:%g] length [%i]', full(min(x)), full(max(x)),length(x));
  else
    x=x(:);
    minmaxstd = sprintf('[%g:%g] size [%s]', full(min(x)), full(max(x)),num2str(size(x)));
  end
  if index==0
    if not(all(m==1 | m==0))
      minmaxstd=[ minmaxstd sprintf(' (per monitor=%g)', mean(m(:))) ];
    end
    minmaxstd=[ minmaxstd sprintf(' sum=%g', sum(iData_private_cleannaninf(x))) ];
  end
  if prod(size(a)) < 1e4
    try
      [s, f] = std(a, index);
      minmaxstd=[ minmaxstd sprintf(' <%g +/- %g>', f,s) ];
    end
  end
  t = sprintf('%6i %15s  %s %s', index, v, l, minmaxstd);
  properties{end+1} = t;
  uimenu(uicm, 'Label', t);
  clear x m
end
% menu About iFit
uimenu(uicm, 'Separator','on','Label', 'About iFit/iData', 'Callback', ...
  [ 'msgbox(''' version(iData,2) sprintf('. Visit <http://ifit.mccode.org>') ''',''About iFit'',''help'')' ]);
% attach contexual menu to plot with UserData storage
ud.properties=properties;
ud.xlabel = xlab;
ud.ylabel = ylab;
ud.zlabel = zlab;
if iscell(titl), titl=titl{1}; end
ud.title  = titl;
ud.name   = char(a);
ud.commands = commandhistory(a);
ud.handle = h;

set(uicm,'UserData', ud);
set(h,   'UIContextMenu', uicm); 

% contextual menu for the axis frame
if ~isempty(get(gca, 'UserData'))
  ud = get(gca, 'UserData');
end
uicm = uicontextmenu;
% menu Duplicate (axis frame/window)
uimenu(uicm, 'Label', 'Duplicate View...', 'Callback', ...
   [ 'tmp_cb.g=gca; tmp_cb.ud=get(gca,''UserData'');' ...
     'tmp_cb.f=figure; tmp_cb.c=copyobj(tmp_cb.g,gcf); ' ...
     'set(tmp_cb.c,''position'',[ 0.1 0.1 0.85 0.8]);' ...
     'set(gcf,''Name'',''Copy of ' char(a) '''); ' ...
     'set(gca,''XTickLabelMode'',''auto'',''XTickMode'',''auto'');' ...
     'set(gca,''YTickLabelMode'',''auto'',''YTickMode'',''auto'');' ...
     'set(gca,''ZTickLabelMode'',''auto'',''ZTickMode'',''auto'');' ...
     'title(tmp_cb.ud.title);', ...
     'xlabel(tmp_cb.ud.xlabel);ylabel(tmp_cb.ud.ylabel); clear tmp_cb;']);
     
if ndims(a) == 1 && ~isfield(ud,'contextual_1d')
  ud.contextual_1d = 1;
end
% menu Toggle all error bars (axis)
if isfield(ud,'contextual_1d') && ud.contextual_1d==1
  uimenu(uicm, 'Label','Toggle All Error Bars', 'Callback', [ ... 
    'tmp_hg = findobj(gca,''type'',''hggroup'');tmp_v=[];'...
    'for tmp_i=1:length(tmp_hg)'...
    'tmp_h=get(tmp_hg(tmp_i),''children'');'...
    'if isempty(tmp_v) ' ...
    'if strcmp(get(tmp_h(2),''visible''),''off''), tmp_v=''on''; else tmp_v=''off''; end;' ...
    'end; set(tmp_h(2),''visible'',tmp_v); clear tmp_h;' ...
    'end; clear tmp_hg tmp_i tmp_v;' ...
   ]);
end
uimenu(uicm, 'Label','Toggle grid', 'Callback','grid');
if ndims(a) >= 2 && ~isfield(ud,'contextual_2d')
  ud.contextual_2d = 1;
end
if isfield(ud,'contextual_2d') && ud.contextual_2d==1
  uimenu(uicm, 'Label','Reset Flat/3D View', 'Callback', [ ...
    '[tmp_a,tmp_e]=view; if (tmp_a==0 & tmp_e==90) view(3); else view(2); end;' ...
    'clear tmp_a tmp_e; lighting none;alpha(1);shading flat;rotate3d off;axis tight;' ]);
  uimenu(uicm, 'Label','Smooth View','Callback', 'shading interp;');
  uimenu(uicm, 'Label','Add Light','Callback', 'light;lighting phong;');
  uimenu(uicm, 'Label','Transparency','Callback', 'alpha(0.7);');
  uimenu(uicm, 'Label','Linear/Log scale','Callback', 'if strcmp(get(gca,''zscale''),''linear'')  set(gca,''zscale'',''log''); else set(gca,''zscale'',''linear''); end');
  uimenu(uicm, 'Label','Toggle Perspective','Callback', 'if strcmp(get(gca,''Projection''),''orthographic'')  set(gca,''Projection'',''perspective''); else set(gca,''Projection'',''orthographic''); end');
else
  uimenu(uicm, 'Label','Reset View', 'Callback','view(2);lighting none;alpha(1);shading flat;axis tight;rotate3d off;');
  uimenu(uicm, 'Label','Linear/Log scale','Callback', 'if strcmp(get(gca,''yscale''),''linear'')  set(gca,''yscale'',''log''); else set(gca,''yscale'',''linear''); end');
end

uimenu(uicm, 'Separator','on','Label', 'About iFit/iData', ...
  'Callback',[ 'msgbox(''' version(iData,2) sprintf('. Visit <http://ifit.mccode.org>') ''',''About iFit'',''help'')' ]);
set(gca, 'UIContextMenu', uicm);
set(gca, 'UserData', ud);

% add rotate/pan/zoom tools to the figure in case java machine is not started
if ~usejava('jvm')
  uicmf = uicontextmenu;
  uimenu(uicmf, 'Label','Zoom on/off', 'Callback','zoom');
  uimenu(uicmf, 'Label','Pan on/off',  'Callback','pan');
  if ndims(a) >= 2
    uimenu(uicmf, 'Label', 'Rotate on/off', 'Callback','rotate3d');
  end
  uimenu(uicmf, 'Label','Legend on/off', 'Callback','legend(gca, ''toggle'',''Location'',''Best'');');
  uimenu(uicmf, 'Label','Print...', 'Callback','printpreview');
  set(gcf, 'UIContextMenu', uicmf);
  set(gcf, 'KeyPressFcn', @(src,evnt) eval('if lower(evnt.Character)==''r'', lighting none;alpha(1);shading flat;axis tight;rotate3d off; zoom off; pan off; end') );
end

try
  set(h,   'Tag',  [ mfilename '_' a.Tag ]);
end
set(gcf, 'Name', char(a));

% labels
if (strfind(method,'hide_ax'))
  % set(gca,'visible','off'); 
  % set(gca,'XTickLabel',[],'XTick',[]); set(gca,'YTickLabel',[],'YTick',[]); set(gca,'ZTickLabel',[],'ZTick',[])
  xlabel(' '); ylabel(' '); zlabel(' ');
  title(T,'interpreter','none');
else
  if ~isempty(xlab), xlabel(xlab,'interpreter','none'); end
  if ~isempty(ylab), ylabel(ylab,'interpreter','none'); end
  if ndims(a) == 3 & ~isempty(clab)
      if iscell(titl)
          titl = { clab , titl{:} };
      else
          titl = { clab, titl };
      end
  end
  if ~isempty(d)
    titl = [ titl ' ''' d '''' ]; 
  end
  title(textwrap(cellstr(T_char),80),'interpreter','none');
end

if ih == 1, hold on; else hold off; end

end
% ============================================================================


