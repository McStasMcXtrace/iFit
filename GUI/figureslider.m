function f=figureslider(varargin)
% FIGURESLIDER cretes by itself a new figure with horizontal and vertical sliders
% and returns its handle. Sliders are hidden if all objects are shown in the figure.
% When the figure gets smaller than the enclosed object set, sliders may be used to 
% shift the figure view.
%
% FIGURESLIDER(H) adds sliders to an existing window.
%
% FIGURESLIDER works the same as FIGURE. You may use GET(H) to obtain properties 
% of the figure.
%
% Example: figureslider('Name','Figure with sliders', 'ToolBar','figure')
%
% See also: figure

  % create figure if none specified
  if isempty(varargin)
    hfig.handle=figure;
  else
    if nargin == 2 && strcmp(varargin{2}, 'resize')
      callback_resize(varargin{1}, [], []);
      return;
    else
      if isnumeric(varargin{1})
        hfig.handle=figure(varargin{1});
        set(hfig.handle,varargin{2:end});
      else
        hfig.handle=figure(varargin{:}); 
      end
    end
  end

  % check if there are already sliders and get their handles
  hfig.slider_h = findobj(hfig.handle,'Tag','Slider_h');
  hfig.slider_v = findobj(hfig.handle,'Tag','Slider_v');
  hfig.center   = findobj(hfig.handle,'Tag','Center');
  hfig.units    = get(hfig.handle, 'Units');
  hfig.sliderwidth=20;

  set(hfig.handle, 'Units','pixels');
  hfig.position  =get(hfig.handle,'Position');
  % horizontal slider
  if isempty(hfig.slider_h)
    hfig.slider_h = uicontrol(hfig.handle,'Style','slider',...
                 'Tag','Slider_h', ...
                 'Position', [1 1 ...
                    hfig.position(3)-hfig.sliderwidth hfig.sliderwidth], ...
                 'Max',100,'Min',0,'Value',50,...
                 'SliderStep',[0.05 0.2],...
                 'UserData',50, ...
                 'Callback',@callback_slider_h);
  end

  % vertical slider
  if isempty(hfig.slider_v)
    hfig.slider_v = uicontrol(hfig.handle,'Style','slider',...
                 'Tag','Slider_v', ...
                 'Position', ... 
                 [hfig.position(3)-hfig.sliderwidth hfig.sliderwidth ...
                   hfig.sliderwidth hfig.position(4)-hfig.sliderwidth], ...
                 'Max',100,'Min',0,'Value',50,...
                 'SliderStep',[0.05 0.2],...
                 'UserData',50, ...
                 'Callback',@callback_slider_v);
  end
  % add a button on lower right size, in between sliders for centering
  if isempty(hfig.center)
    hfig.center = uicontrol(hfig.handle,'Style','pushbutton',...
                 'Tag','Center', ...
                 'Position', ... 
                 [hfig.position(3)-hfig.sliderwidth 1 ...
                   hfig.sliderwidth hfig.sliderwidth], ....
                 'ToolTip','Click here to center figure contents', ...
                 'Callback',@callback_center);
  end

  % set resize callback so that sliders may be hidden/shown 
  % depending on size of inner objects to show
  resizeFcn = get(hfig.handle, 'ResizeFcn');
  if ~isempty(resizeFcn)
    if isempty(findstr(char(resizeFcn), 'figureslider(gcf, ''resize'');'))
      set(hfig.handle, 'ResizeFcn', [  'figureslider(gcf, ''resize''); ; ' char(resizeFcn) ';' ]);
    end
  else
    set(hfig.handle, 'ResizeFcn', 'figureslider(gcf, ''resize'');');
  end

  callback_resize(hfig.handle);  % setup sliders

  % restore initial figure Unit settings
  set(hfig.handle, 'Units', hfig.units);
  
  f=hfig.handle;

% ========== private functions ================================================
function callback_center(hObject, eventdata, handles)
% center figure on enclosing objects
  hfig.handle   = gcf;
  hfig.slider_h = findobj(hfig.handle,'Tag','Slider_h');
  hfig.slider_v = findobj(hfig.handle,'Tag','Slider_v');
  set(hfig.slider_h, 'Value',50);
  callback_slider_h(hfig.slider_h);
  set(hfig.slider_v, 'Value',50);
  callback_slider_v(hfig.slider_v);

function callback_slider_h(hObject, eventdata, handles)
% move enclosing objects if not in units=nomalized 
  hfig.handle   = gcf;
  hfig.children = get(hfig.handle, 'Children');
  hfig.units    = get(hfig.handle, 'Units');
  hfig.slider_h = findobj(hfig.handle,'Tag','Slider_h');
  hfig.slider_v = findobj(hfig.handle,'Tag','Slider_v');
  hfig.center   = findobj(hfig.handle,'Tag','Center');
  set(hfig.handle, 'Units','pixels');
  extension = callback_getchildrenextension(hfig.handle);
  
  for index=1:length(hfig.children)
    this.handle = hfig.children(index);
    this.pos    = get(hfig.children(index), 'Position');
    try
        this.units  = get(hfig.children(index), 'Units');
    catch
        this.units  = [];
    end
    % skip sliders
    if isempty(this.units) | this.handle == hfig.slider_h | this.handle == hfig.slider_v | this.handle == hfig.center
      continue; 
    end
    if strcmp(this.units,'normalized'), continue; 
    else set(hfig.children(index), 'Units', 'pixels');
    end
    step = [ extension(3)*(get(hObject,'Value')-get(hObject,'UserData'))/100 0 0 0 ];
    set(hfig.children(index), 'Position', this.pos-step);
    set(hfig.children(index), 'Units', this.units);
  end
  set(hObject,'UserData', get(hObject,'Value'));
  
  set(hfig.handle, 'Units', hfig.units);

function callback_slider_v(hObject, eventdata, handles)
% move enclosing objects if not in units=nomalized 
  hfig.handle   = gcf;
  hfig.children = get(hfig.handle, 'Children');
  hfig.units    = get(hfig.handle, 'Units');
  hfig.slider_h = findobj(hfig.handle,'Tag','Slider_h');
  hfig.slider_v = findobj(hfig.handle,'Tag','Slider_v');
  hfig.center   = findobj(hfig.handle,'Tag','Center');
  set(hfig.handle, 'Units','pixels');
  extension = callback_getchildrenextension(hfig.handle);
  
  for index=1:length(hfig.children)
    this.handle = hfig.children(index);
    this.pos    = get(hfig.children(index), 'Position');
    try
        this.units  = get(hfig.children(index), 'Units');
    catch
        this.units  = [];
    end
    % skip sliders
    if isempty(this.units) | this.handle == hfig.slider_h | this.handle == hfig.slider_v | this.handle == hfig.center
      continue; 
    end
    if strcmp(this.units,'normalized'), continue; 
    else set(hfig.children(index), 'Units', 'pixels');
    end
    step = [ 0 extension(4)*(get(hObject,'Value')-get(hObject,'UserData'))/100 0 0 ];
    set(hfig.children(index), 'Position', this.pos-step);
    set(hfig.children(index), 'Units', this.units);
  end
  set(hObject,'UserData', get(hObject,'Value'));
  
  set(hfig.handle, 'Units', hfig.units);

function callback_resize(hObject, eventdata, handles)
  hfig.handle   = hObject;
  hfig.units    = get(hfig.handle, 'Units');
  set(hfig.handle, 'Units','pixels');
  hfig.position = get(hfig.handle,'Position');
  hfig.slider_h = findobj(hfig.handle,'Tag','Slider_h');
  hfig.slider_v = findobj(hfig.handle,'Tag','Slider_v');
  hfig.center   = findobj(hfig.handle,'Tag','Center');

  hfig.sliderwidth=20;

  % move sliders at figure edges
  set(hfig.slider_h, 'Position', ...
    [1 1 ...
     hfig.position(3)-hfig.sliderwidth hfig.sliderwidth]);
  set(hfig.slider_v, 'Position', ... 
    [hfig.position(3)-hfig.sliderwidth hfig.sliderwidth ...
     hfig.sliderwidth hfig.position(4)-hfig.sliderwidth]);
  set(hfig.center, 'Position', ... 
    [hfig.position(3)-hfig.sliderwidth 1 ...
     hfig.sliderwidth hfig.sliderwidth]);

  % hide or show sliders, depending on size of uipanel vs. figure
  extension = callback_getchildrenextension(hfig.handle);

  if extension(3) <= hfig.position(3)
    set(hfig.slider_h, 'Visible','off');
  else
    set(hfig.slider_h, 'Visible','on');
  end

  if extension(4) <= hfig.position(4)
    set(hfig.slider_v, 'Visible','off');
  else
    set(hfig.slider_v, 'Visible','on');
  end
  
  if extension(3) <= hfig.position(3) & extension(4) <= hfig.position(4)
    set(hfig.center, 'Visible','off');
  else
    set(hfig.center, 'Visible','on');
  end

  set(hfig.handle, 'Units', hfig.units);

function extension=callback_getchildrenextension(hObject, eventdata, handles)
% compute the total extension on enclosed objects
  hfig.handle   = hObject;
  hfig.units    = get(hfig.handle, 'Units');
  set(hfig.handle, 'Units','pixels');
  hfig.children = get(hfig.handle, 'Children');
  hfig.slider_h = findobj(hfig.handle,'Tag','Slider_h');
  hfig.slider_v = findobj(hfig.handle,'Tag','Slider_v');
  hfig.center   = findobj(hfig.handle,'Tag','Center');

  min_x=Inf; d_x=-Inf; min_y=Inf; d_y=-Inf;
  for index=1:length(hfig.children)
    this.handle = hfig.children(index);
    this.pos    = get(hfig.children(index), 'Position');
    try
        this.units  = get(hfig.children(index), 'Units');
    catch
        this.units  = [];
    end
    % skip sliders
    if isempty(this.units) | this.handle == hfig.slider_h | this.handle == hfig.slider_v | this.handle == hfig.center
      continue; 
    end
    if strcmp(this.units, 'normalized'), continue;
    else set(hfig.children(index), 'Units', 'pixels');
    end
    % get extension of children
    min_x = min(this.pos(1), min_x);        % origin of object
    min_y = min(this.pos(2), min_y);
    d_x   = max(this.pos(1)+this.pos(3), d_x);  % maximum extension of object
    d_y   = max(this.pos(2)+this.pos(4), d_y);
    set(hfig.children(index), 'Units', this.units);
  end

  extension = [ min_x min_y d_x-min_x d_y-min_y ];

  set(hfig.handle, 'Units', hfig.units);

