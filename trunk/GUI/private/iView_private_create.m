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
    iview(instance, 'config_save');  % first instance: we save the configuration
  end
  instance_list = [ instance_list(:) ; instance ];
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
    instance = figureslider('ResizeFcn', 'iview(gcf, ''instance_resize'');', ...
    'CloseRequestFcn','iview(gcf, ''instance_close'');');  % new figure with sliders
  else
    instance = figureslider(instance,'ResizeFcn', 'iview(gcf, ''instance_resize'');', ...
    'CloseRequestFcn','iview(gcf, ''instance_close'');');
  end
  set(instance, 'MenuBar','none', 'ToolBar','none');
  set(instance, 'PaperPositionMode','auto');
  set(instance, 'Tag','iView_instance');
  set(instance, 'Name',[ 'iView #' num2str(instance) ]);
  set(instance, 'HandleVisibility', 'callback','Interruptible','off');
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
  if isfield(config,'PaperType')
    set(instance, 'PaperType',       config.PaperType);
  else
    config.PaperType = get(instance, 'PaperUnits');
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
  uimenu(file, 'Label', '&New window','Callback','iview(''instance_new'');','Accelerator','n');
  uimenu(file, 'Label', 'New data set', 'Callback','iview(''data_new'');');
  uimenu(file, 'Label', '&Open...','Callback','iview(gcf, ''data_load'');','Accelerator','o');
  uimenu(file, 'Label', '&Save', 'Callback', 'iview(gcf, ''data_save'',''selection'');','Accelerator','s', 'Separator','on');
  uimenu(file, 'Label', 'Save as...', 'Callback', 'iview(gcf, ''data_saveas'',''selection'');','Accelerator','s');
  uimenu(file, 'Label', 'Save configuration', 'Callback', 'iview(gcf, ''config_save'');');
  uimenu(file, 'Label', 'Page setup...', 'Callback', 'pagesetupdlg(gcf);', 'Separator','on');
  uimenu(file, 'Label', '&Print...',      'Callback', 'printdlg(gcf);','Accelerator','p');
  uimenu(file, 'Label', 'Preferences...', 'Enable','off');
  uimenu(file, 'Label', 'Close &window', 'Callback', 'iview(gcf, ''instance_close'');', 'Separator','on','Accelerator','w');
  uimenu(file, 'Label', '&Exit', 'Callback', 'iview(gcf, ''exit'');','Accelerator','q'); % quit application
  
  edit = uimenu(instance, 'Label', '&Edit');
  uimenu(edit, 'Label', 'Cut', 'Enable','off','Accelerator','x');
  uimenu(edit, 'Label', 'Copy', 'Enable','off','Accelerator','c');
  uimenu(edit, 'Label', 'Paste', 'Enable','off','Accelerator','v');
  uimenu(edit, 'Label', 'Select &all', 'Callback', 'iview(gcf, ''data_select_all'');', 'Separator','on','Accelerator','a');
  uimenu(edit, 'Label', '&Deselect all', 'Callback', 'iview(gcf, ''data_deselect_all'');','Accelerator','d');
  uimenu(edit, 'Label', '&Find...', 'Enable','off', 'Separator','on','Accelerator','f'); % dialog to find match, and select result
  uimenu(edit, 'Label', '&Rename window...', 'Callback', 'iview(gcf, ''instance_rename'');');

  % create dynamic menu (from config)

  % create static contextual menu
  cmenu = uicontextmenu('Parent',instance);
  set(instance, 'UIContextMenu', cmenu);
  uimenu(cmenu, 'Label', 'New window','Callback','iview(''instance_new'');');
  uimenu(cmenu, 'Label', 'New data set', 'Callback','iview(''data_new'');');
  uimenu(cmenu, 'Label', 'Open...', 'Callback','iview(gcf, ''data_load'');');
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
  uimenu(documents, 'Label', '&Open data (plot)', 'Tag','Static', 'Callback', 'iview(gcf, ''data_open'', ''selection'');');
  uimenu(documents, 'Label', '&Edit data...', 'Tag','Static', 'Enable','off');  % edit content/properties/axes/signal/alias...
  uimenu(documents, 'Label', '&Properties...', 'Tag','Static', 'Callback', 'iview(gcf, ''data_properties'', ''selection'');','Accelerator','i');
  sortby = uimenu(documents, 'Label', 'Sort by', 'Tag','Static'); % Size, Date, Title, Label, Tag
  	uimenu(sortby, 'Label', 'Title', 'Tag','Static', 'Callback','iview(gcf, ''instance_sort'',''Title'');'); 
  	uimenu(sortby, 'Label', 'Date', 'Tag','Static', 'Callback','iview(gcf, ''instance_sort'',''Date'');'); 
  	uimenu(sortby, 'Label', 'Size', 'Tag','Static', 'Callback','iview(gcf, ''instance_sort'',''Size'');'); 
  	uimenu(sortby, 'Label', 'Label', 'Tag','Static', 'Callback','iview(gcf, ''instance_sort'',''Label'');'); 
  	uimenu(sortby, 'Label', 'Tag (unique ID)', 'Tag','Static', 'Callback','iview(gcf, ''instance_sort'',''Tag'');');
  uimenu(documents, 'Label', 'Rename...', 'Tag','Static', 'callback',	'iview(gcf, ''data_rename'', ''selection'');');
  uimenu(documents, 'Label', 'Set Label as', 'Tag','Static', 'Separator','on');
  uimenu(documents, 'Label', 'Define new label...', 'Tag','Static', 'Callback','iview(gcf, ''data_label'', ''add'');');
  uimenu(documents, 'Label', 'Remove label...', 'Tag','Static', 'Callback', 'iview(gcf, ''data_label'', ''remove'');');
  uimenu(documents, 'Separator','on', 'Label', 'Delete selection', 'Tag','Static', 'Callback', 'iview(gcf, ''data_close'', ''selection'');');
  
  help=uimenu(instance, 'Label', 'Help');
  uimenu(help, 'Label', 'Contents', 'Enable','off');
  uimenu(help, 'Label', 'Contacts', 'Callback', 'iview(gcf, ''contacts'');');
  uimenu(help, 'Label', 'About iView', 'Callback', 'iview(gcf, ''about'');', 'Separator','on');
  
  % install mouse/keyboard event handler
  set(instance, 'ButtonDownFcn', 'iview(gcf, ''mouse_down'', gcf);');
  % keyboard is handled through menu accelerators
  
  movegui(instance); % make sure the window is visible on screen
  
  % create Label items and update Documents menu
  iView_private_documents(instance);
  
  config.Position = get(instance, 'Position');

function config=iView_private_create_storage(instance, config)
% Create storage area (AppData) for iData and Config
  setappdata(instance, 'Data', []);
  
