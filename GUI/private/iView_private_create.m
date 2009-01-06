function [instance,config]=iView_private_create(instance, config);
% iView_private_create Create or raise an iView instance

if nargin == 0, instance=[]; end
if nargin == 1, config = []; end

% find all existing iView instances
instance_list=findall(0,'Tag','iView_instance');

% check instance
create_new_instance  = 0;

if isempty(instance_list)
  create_new_instance  = 1;   % no instance exists: create one
elseif isempty(instance)
  instance=instance_list(1);  % use first existing instance
elseif isempty(find(instance == instance_list))
  create_new_instance  = 1;   % invalid specified instance: create one
end

if isempty(config), config=iView_private_config(instance, 'load'); end

% create interface and associated storage areas (iData, config)
if create_new_instance
  % create figure and Name
  [instance, config]=iView_private_create_interface(instance, config);
  % create storage areas
  iView_private_create_storage(instance, config);
  if ~length(instance_list)
    disp([ '% ' datestr(now) ' Starting iView' ]);
    iView_private_config(0, 'save', config);
  end
  instance_list = [ instance_list instance ];
end

% instance is valid. Raise it, update figure Names and return.
% a '*' is prepended to figure Name to indicate active iView instance
for index=1:length(instance_list)
  this_instance= instance_list(index);
  this_name    = get(this_instance, 'Name');
  
  if this_name(1) ~= '*' & this_instance == instance
      this_name= [ '*' this_name ]; % add '*' to Name
  elseif this_name(1) == '*' & this_instance ~= instance
      this_name(1) = '';            % remove '*'
  end
  set(this_instance, 'Name', this_name);
end % for

% raise figure and return
figure(instance);

% =============================================================================
%                    INLINE functions
% =============================================================================
% iView_private_create_interface
% iView_private_create_storage

function [instance, config]=iView_private_create_interface(instance, config)
  % Build main iView interface (window, menus, ...)

  % create main windows with sliders
  if isempty(instance)
    instance = figureslider('ResizeFcn', 'iview(gcf, ''resize'');');  % new figure with sliders
  else
    instance = figureslider(instance,'ResizeFcn', 'iview(gcf, ''resize'');');
  end
  set(instance, 'MenuBar','none', 'ToolBar','none');
  set(instance, 'Tag','iView_instance');
  set(instance, 'Name',[ 'iView #' num2str(instance) ]);
  set(instance, 'HandleVisibility', 'callback');
  set(instance, 'NextPlot','new', 'NumberTitle','on', 'BackingStore','on');
  
  % install event handler for mouse
  %set(instance, 'WindowButtonDownFcn',...
	%					'iview(gcf, ''mouse_down'', gco);');
  
  % remove horizontal slider
  slider_h = findobj(instance,'Tag','Slider_h');
  delete(slider_h);
  
  % customize interface from config
  if isfield(config,'Color')
    set(instance, 'Color',            config.Color);
  else
    config.Color = get(instance, 'Color');
  end
  if isfield(config,'FileName')
    set(instance, 'FileName',         config.FileName);
  else
    config.FileName = get(instance, 'FileName');
  end
  if isfield(config,'PaperUnits')
    set(instance, 'PaperUnits',       config.PaperUnits);
  else
    config.PaperUnits = get(instance, 'PaperUnits');
  end
  if isfield(config,'PaperOrientation')
    set(instance, 'PaperOrientation', config.PaperOrientation);
  else
    config.PaperOrientation = get(instance, 'PaperOrientation');
  end
  if isfield(config,'PaperPosition')
    set(instance, 'PaperPosition', config.PaperPosition);
  else
    config.PaperPosition = get(instance, 'PaperPosition');
  end
  
  if isfield(config,'Size')
    pos      = get(instance, 'Position');
    pos(3:4) = config.Size;
    set(instance, 'Position', pos);
  end
  if isfield(config,'Position')
    set(instance, 'Position', config.Position);
  else
    config.Position = get(instance, 'Position');
  end

  % create static menus (based on default Figure menus)
  file = uimenu(instance, 'Label', '&File');
  uimenu(file, 'Label', '&New window','Callback','iview(''new'');','Accelerator','n');
  uimenu(file, 'Label', 'New data set', 'Callback','iview(''new_data'');');
  uimenu(file, 'Label', '&Open...','Callback','iview(gcf, ''load'');','Accelerator','o');
  uimenu(file, 'Label', '&Save', 'Callback', 'iview(gcf, ''save_data'');','Accelerator','s', 'Separator','on');
  uimenu(file, 'Label', 'Save as...', 'Callback', 'iview(gcf, ''saveas_data'');','Accelerator','s');
  uimenu(file, 'Label', 'Save configuration', 'Callback', 'iview(gcf, ''save_config'');');
  uimenu(file, 'Label', 'Page setup...', 'Callback', 'pagesetupdlg(gcf);', 'Separator','on');
  uimenu(file, 'Label', '&Print...',      'Callback', 'printdlg(gcf);','Accelerator','p');
  uimenu(file, 'Label', 'Preferences...', 'Enable','off');
  uimenu(file, 'Label', 'Close &window', 'Callback', 'iview(gcf, ''close'');', 'Separator','on','Accelerator','w');
  uimenu(file, 'Label', '&Exit', 'Callback', 'iview(gcf, ''exit'');','Accelerator','q'); % quit application
  
  edit = uimenu(instance, 'Label', '&Edit');
  uimenu(edit, 'Label', 'Cut', 'Enable','off','Accelerator','x');
  uimenu(edit, 'Label', 'Copy', 'Enable','off','Accelerator','c');
  uimenu(edit, 'Label', 'Paste', 'Enable','off','Accelerator','v');
  uimenu(edit, 'Label', 'Select &all', 'Callback', 'iview(gcf, ''select_all'');', 'Separator','on','Accelerator','a');
  uimenu(edit, 'Label', '&Deselect all', 'Callback', 'iview(gcf, ''deselect_all'');','Accelerator','d');
  uimenu(edit, 'Label', '&Find...', 'Enable','off', 'Separator','on','Accelerator','f'); % dialog to find match, and optionally select result
  sort = uimenu(edit, 'Label', '&Sort as...'); % Date, Size, Name, Label
  
  view = uimenu(instance, 'Label', '&View', 'Enable','off'); 
  %uimenu(view, 'Label', 'Menu'); % only in uicontext menu
  uimenu(view, 'Label', 'Refresh');
  uimenu(view, 'Label', 'Toolbar', 'Separator','on');
  uimenu(view, 'Label', 'Icons'); % toggle Icon view
  uimenu(view, 'Label', 'Icon size...');

  % create dynamic menu (from config)

  % create static contextual menu
  cmenu = uicontextmenu('Parent',instance);
  set(instance, 'UIContextMenu', cmenu);
  uimenu(cmenu, 'Label', 'New data set', 'Callback','iview(''new_data'');');
  uimenu(cmenu, 'Label', 'Open...', 'Callback','iview(gcf, ''load'');');
  % uicontext menu on background:
  %   new window
  %   load
  %   save 
  %   select all
  %   deselect all
  %   properties
  %   paste
  %   align

  % create dynamic contextual menu (from config)
  
  % menus that must be at the right side
  documents=uimenu(instance, 'Label', 'Documents','Tag','Documents');
  uimenu(documents, 'Label', '&Edit data...', 'Tag','Static', 'Enable','off');
  uimenu(documents, 'Label', '&Properties', 'Tag','Static', 'Callback', 'iview(gcf, ''properties'', ''selection'');','Accelerator','i');
  uimenu(documents, 'Label', 'Delete', 'Tag','Static', 'Callback', 'iview(gcf, ''close_data'', ''selection'');');
  uimenu(documents, 'Label', 'Move to new window', 'Tag','Static', 'Enable','off');
  %uimenu(documents, 'Separator','on');
  
  help=uimenu(instance, 'Label', 'Help');
  uimenu(help, 'Label', 'Contents', 'Enable','off');
  uimenu(help, 'Label', 'Contacts', 'Enable','off');
  uimenu(help, 'Label', 'About iView', 'Callback', 'iview(gcf, ''about'');', 'Separator','on');
  
  % install mouse/keyboard event handler
  set(instance, 'ButtonDownFcn', 'iview(gcf, ''mouse_down'', gcf);');
  % keyboard is handled throught menu accelerators
  
  movegui(instance); % make sure the window is visible on screen
  
  config.Position = get(instance, 'Position');

function config=iView_private_create_storage(instance, config)
% Create storage area (AppData) for iData and Config
  setappdata(instance, 'Data', []);
  
