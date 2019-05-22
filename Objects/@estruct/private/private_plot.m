function h =private_plot(a, method, varargin)

h      = [];
if isempty(a)
  return; 
end

if isempty(method)
  if isvector(a) > 3
    method='scatter';
  else
    method='plot';
  end
end

% clean method string from the plot type and supported options not to be passed to matlab plot commands
if ischar(method)
  toremove='plot3 stem3 scatter3 scatter stem plot mesh surf waterfall tight auto hide view2 view3 transparent axis hide_err hide_errorbars hide_error contour contour3 surfc surfl contourf pcolor median mean half slice flat interp faceted light clabel colorbar shifted hide_axes painters zbuffer whole full legend replace update grid';
  toremove=strread(toremove,'%s','delimiter',' ');
  this_method = method;
  for index=1:length(toremove)
    this_method= regexprep(this_method,[ '\<' toremove{index} '\>' ], '');
  end
else
  this_method = method;
end
this_method = strtrim(this_method);

% plot an array of objects =====================================================
if numel(a) > 1
  h = private_plot_array(a, method, this_method, varargin{:});
  return
end % plot array

% plot a single object
method = lower(method);

% check if the object is not too large, else rebin accordingly
if prod(size(a)) > 1e6 
  if isempty([ strfind(method,'whole') strfind(method,'full') ])
    warning([ mfilename ': Object ' a.Tag ' "' a.Name '" is large (numel=' num2str(prod(size(a))) ...
      ').\n\tNow rebinning for display purposes with e.g. a=reducevolume(a);' ...
      '\n\tUse e.g plot(a, ''whole'') to plot the whole data set and be able to zoom tiny regions.' ]);
    tag = a.Tag;
    a=reducevolume(a);
    a.Tag = tag;
  else
    method = [ method ' opengl' ];
  end
end
zlab = '';

% replace/update existing plot
if ~isempty(strfind(method,'update')) || ~isempty(strfind(method,'replace'))
  h = [ findall(0, 'Tag', [ mfilename '_' a.Tag ]) ...
        findall(0, 'Tag', [ mfilename '_contextmenu_' a.Tag ]) ];
  if ~isempty(h)
    % search parent figure
    ax = get(h(1),'Parent');
    try
        while ~strcmp(get(ax,'Type'),'figure')
            ax = get(ax,'Parent');
        end
    end
    figure(ax);
    delete(h);
  end
elseif ~isempty(strfind(method,'figure'))
  figure;
end

% possibly select Rendered prior to start plotting
if exist('feature') && ~feature('ShowFigureWindows')
  set(gcf,'Renderer','painters')
elseif ~isempty(strfind(method,'opengl'))   % faster for large data sets
	set(gcf,'Renderer','OpenGL')
elseif ~isempty(strfind(method,'painters'))
	set(gcf,'Renderer','painters')
elseif ~isempty(strfind(method,'zbuffer'))
	set(gcf,'Renderer','zbuffer');
elseif ismac
  set(gcf,'Renderer','painters'); % default for MacOS which do not support OpenGL
end

% ==============================================================================
ret = 0;

m = []; mp = []; mv = []; names = []; name = [];
% get Model,etc... when found in the Dataset

if isfield(a, 'Model')
  m = get(a, 'Model');
elseif ~isempty(findfield(a, 'Model'))
  m = get(a, findfield(a, 'Model', 'cache first'));
end

if isa(m, 'iFunc') && ~isempty(m)
  % get the parameter values as a struct
  mp    = m.ParameterValues;
  names = m.Parameters; names = names(:);
  name  = m.Name;
end


if isfield(a, 'ModelValue')
  mv = get(a, 'ModelValue');
  if ~strcmp(getalias(a,'Signal'), 'ModelValue')
    set(a, 'ModelValue', []);  % avoid recursive loop
  end
elseif ~isempty(findfield(a, 'ModelValue'))
  mv = get(a, findfield(a, 'ModelValue', 'cache first'));
end

if isempty(mp)
  if isfield(a, 'ModelParameters')
    mp = get(a, 'ModelParameters');
  elseif isfield(mv, 'ModelParameters')
    mp = get(mv, 'ModelParameters');
  end
end
% make it a structure
if ~isempty(mp) && numel(names) == numel(mp)
  mp = cell2struct(num2cell(mp(:)),strtok(names(:)));
end

if ~isempty(mv) && isa(mv, class(a)) ...
  && isfield(mv, 'ModelValue') && ~strcmp(getalias(mv,'Signal'), 'ModelValue')
  set(mv,'ModelValue', []);
end

switch ndims(a) % handle different plotting methods depending on the object dimensionality
case 0
  h=[]; 
  return;
case 1  % vector type data (1 axis + signal) -> plot
  [h, xlab, ylab, ret] = private_plot_1d(a, method, this_method, varargin{:}); % in private
case 2  % surface type data (2 axes+signal) -> surf or plot3
  [h, xlab, ylab, zlab] = private_plot_2d(a, method, this_method, varargin{:}); % in private
otherwise % 3d data sets: volumes
  if ndims(a) > 3
    % reduce dimensions: look for the axes which have the largest extension
    extent = 1:ndims(a);
    for index=1:ndims(a)
      this_axis = getaxis(a, index);
      extent(index) = max(this_axis(:)) - min(this_axis(:));
    end
    % now get the largest extents
    [extent,extent_index] = sort(extent);
    extent_index = extent_index(end:-1:1);  % descending order
    extent       = extent(end:-1:1);
    % we use extend_index(1:3)
    % we keep axes [1:2 ndims(a)] (this allows e.g. to plot S(q,w))
    if isvector(a) > 1 % event data set
      to_remove = [];
      for index=2:ndims(a)
        if isempty([ strfind(method,'whole') strfind(method,'full') ])
          if index <= 4 && extent(index) > extent(1)*1e-2, continue; end
        elseif index <= 3, continue; end
        to_remove = [ to_remove extent_index(index) ];
      end
      a = rmaxis(a,to_remove);
    else
      sz = size(a); sz(extent_index(4:end)) = 1;
      warning([ mfilename ': Reducing ' num2str(ndims(a)) '-th dimensional data ' a.Tag ' "' a.Name '" to 3D with a=resize(a, ' mat2str(sz) ')' ]);
      a = squeeze(resize(a, sz));
    end
  end
  [h, xlab, ylab, zlab, ret] = private_plot_3d(a, method, this_method, varargin{:}); % in private
end % switch

if ret
  return
end

% tune the rendering of the plot ===============================================
if ~isempty(strfind(method,'flat'))
  shading flat
elseif ~isempty(strfind(method,'interp'))
  shading interp
elseif ~isempty(strfind(method,'faceted'))
  shading faceted
end
if ~isempty(strfind(method,'transparent')) || ~isempty(strfind(method,'alpha'))
  alpha(0.7);
end
if ~isempty(strfind(method,'light')) && isempty(findobj(gca,'Type','light'))
  light;
end
if ~isempty(strfind(method,'view2'))
  view(2);
end
if ~isempty(strfind(method,'view3'))
  view(3);
end
if ~isempty(strfind(method,'tight'))
  axis tight
end
if ~isempty(strfind(method,'auto'))
  axis auto
end
if ~isempty(strfind(method,'colorbar'))
  cb = colorbar;
  title(cb, label(a, 0));
end

% add a UIcontextMenu so that right-click gives info about the object plot
T   = a.Name;
if isempty(T), T=label(a,0); end
if ~ischar(T), T=char(T); end
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
  if strcmp(a.Label, a.DisplayName)
      if ~isempty(title(a)), a.DisplayName=title(a);
      else a.DisplayName=fS; end
  end
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
T0 = T; % original title, full.

if length(T) > 23, T=[ T(1:20) '...' ]; end
if length(S)+length(d) < 30,
  d = [ d ' ' T ];
end


% install the contextual menu
private_plot_contextmenu(a, h, xlab, ylab, zlab, T, S, d, cmd, mp, name);

% assign some settings to graphics handles, including Children
children = get(h,'Children');
if ~isempty(strtrim(d)), displayname = d; else displayname=[ T a.Tag ' <' S '>' ]; end
try
  set(h , 'DisplayName', displayname);
  set(children , 'DisplayName', displayname);
end

try
  set(h,          'Tag',  [ mfilename '_' a.Tag ]);
  set(children,   'Tag',  [ mfilename '_' a.Tag ]);
end
set(gcf, 'Name', char(a));

% labels
if ~isempty(strfind(method,'hide_ax'))
  % set(gca,'visible','off'); 
  % set(gca,'XTickLabel',[],'XTick',[]); set(gca,'YTickLabel',[],'YTick',[]); set(gca,'ZTickLabel',[],'ZTick',[])
  xlabel(' '); ylabel(' '); zlabel(' ');
  title(T,'interpreter','none');
else
  if ~isempty(xlab), xlabel(xlab); end
  if ~isempty(ylab), ylabel(ylab); end

  title(textwrap(cellstr(char(T0)),80),'interpreter','none');
end

if (strfind(method,'legend'))
  legend(h);
end

% display model value if available
if isa(mv, class(a)) && ~isempty(mv) && ~all(all(isnan(mv))) && (ndims(mv) <= 2 || isvector(mv) > 1)
  hold on
  axis(axis); % fix plot limits
  h2 = [];
  if isvector(mv) > 1
    h2 = plot(mv,'r--');
  else
    switch ndims(mv)
    case 1
      h2 = plot(mv,'r--');
    case 2
      h2 = contour3(mv);
    end
  end
  h = [ h(:) ;  h2(:) ];
end

% ============================================================================



