function instance=iView_private_create(instance, config);
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
  
  % customize interface from config
  if isfield(config,'Color')
    set(instance, 'Color',            config.Color);
  end
  if isfield(config,'FileName')
    set(instance, 'FileName',         config.FileName);
  end
  if isfield(config,'PaperUnits')
    set(instance, 'PaperUnits',       config.PaperUnits);
  end
  if isfield(config,'PaperOrientation')
    set(instance, 'PaperOrientation', config.PaperOrientation);
  end
  if isfield(config,'PaperPosition')
    set(instance, 'PaperPosition',config.PaperPosition);
  end
  
  if isfield(config,'Position')
    set(instance, 'Position', config.Position); 
  end
  if isfield(config,'Size')
    pos = get(instance, 'Position');
    pos(3:4) = config.Size;
    set(instance, 'Position', pos);
  end

  % create static menus (based on default Figure menus)
  file = uimenu(instance, 'Label', 'File');
  uimenu(file, 'Label', 'New window','Callback','iview(''new'');');
  uimenu(file, 'Label', 'Open...','Callback','iview(gcf, ''load'');');
  uimenu(file, 'Label', 'Save');
  uimenu(file, 'Label', 'Save as...');
  uimenu(file, 'Label', 'Page setup...', 'Callback', 'pagesetupdlg(gcf);', 'Separator','on');
  uimenu(file, 'Label', 'Print...',      'Callback', 'printdlg(gcf);');
  uimenu(file, 'Label', 'Preferences...');
  uimenu(file, 'Label', 'Close', 'Callback', 'iview(gcf, ''close'');', 'Separator','on');
  uimenu(file, 'Label', 'Exit');
  
  edit = uimenu(instance, 'Label', 'Edit');
  uimenu(edit, 'Label', 'Cut');
  uimenu(edit, 'Label', 'Copy');
  uimenu(edit, 'Label', 'Paste');
  uimenu(edit, 'Label', 'Find...', 'Separator','on'); % dialog to find match, and optionally select result
  uimenu(edit, 'Label', 'Select all');
  uimenu(edit, 'Label', 'Deselect all');
  uimenu(edit, 'Label', 'Properties...');
  
  view = uimenu(instance, 'Label', 'View');
  %uimenu(view, 'Label', 'Menu'); % only in uicontext menu
  uimenu(view, 'Label', 'Refresh');
  uimenu(view, 'Label', 'Toolbar', 'Separator','on');
  uimenu(view, 'Label', 'Icons'); % toggle Icon view
  uimenu(view, 'Label', 'Icon size...');

  % create dynamic menu (from config)

  % create static contextual menu
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
  uimenu(documents, 'Label', 'Save all as...');
  uimenu(documents, 'Label', 'Clear all...');
  uimenu(documents, 'Label', 'Move to new window');
  uimenu(documents, 'Separator','on');
  
  help=uimenu(instance, 'Label', 'Help');
  uimenu(help, 'Label', 'Contents');
  uimenu(help, 'Label', 'Contacts');
  uimenu(help, 'Separator','on');
  uimenu(help, 'Label', 'About iView');
  
  movegui(instance); % make sure the window is visible on screen
  
  config.Position = get(instance, 'Position');

function config=iView_private_create_storage(instance, config)
% Create storage area (AppData) for iData and Config
  setappdata(instance, 'Config', config);
  setappdata(instance, 'Data', []);
  
  
