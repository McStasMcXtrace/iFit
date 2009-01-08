function config=iView_private_config(instance, action, config)
% handle iView instance config
% action: load, save

switch action
case 'load'
  config=iView_private_config_load(0);
case 'save'
  if nargin<3, config=[]; end
  if isempty(config), config=iView_private_config_load(0); end
  config=iView_private_config_save(config);
  if ~isappdata(0, 'iView_Config')  % first time we save the config ? (e.g. starting iView)
    disp([ '% Saved iView configuration into ' config.FileName ]);
  end
  setappdata(0, 'iView_Config', config);
otherwise
  disp([' Unknown action ' action ' in ' mfilename ]);
end

% =============================================================================
%                    INLINE functions
% =============================================================================
% iView_private_config_load
% iView_private_config_save
% iView_private_config_save_rec

function config=iView_private_config_load(instance)
% iView_private_config_load load configuration from config file, or existing window, or default
  config=[];
  
  % search config file
  if exist('iView_ini')
    config = iView_ini;
  elseif exist(fullfile(prefdir, 'iView.ini'), 'file')
    % there is an iLoad_ini in the Matlab preferences directory: read it
    file = fullfile(prefdir, 'iView.ini');
    fid = fopen(file, 'r');
    content = fread(fid, Inf, 'uint8=>char');
    fclose(fid);
    % evaluate content of file
    config=[]; eval(content(:)'); % this make a 'config' variable
    config.FileName = file;
  end

  if isempty(config)
    % is there an existing config ? If yes, re-use it
    if ishandle(instance)
      if isappdata(instance, 'iView_Config')
        config=getappdata(instance, 'iView_Config');
      end
    else % iView instance does not exist, search for an other instance to get config
      % find all existing iView instances
      instance_list=findall(0,'Tag','iView_instance');
      if ~isempty(instance_list)
        instance=instance_list(1);  % use first existing instance
        config=getappdata(instance, 'iView_Config');
      end
    end
  end
  
  % make sure some fields exist and apply default if required
  if ~isfield(config,'FileName'),         config.FileName         = 'iView.ini'; end
  if ~isfield(config,'PaperUnits'),       config.PaperUnits       = 'centimeters'; end
  if ~isfield(config,'PaperOrientation'), config.PaperOrientation = 'landscape'; end
  if ~isfield(config,'PaperPosition'),    config.PaperPosition    = [5 5 20 15]; end
  if ~isfield(config,'PaperType'),        config.PaperType        = 'A4'; end
  if ~isfield(config,'Color'),            config.Color            = [1 1 1];  end
  if ~isfield(config,'IconSize'),         config.IconSize         = 64; end % in pixels. Recommended: 64
  if ~isfield(config,'IconWidth'),        config.IconWidth        = 2; end  % real icon width = IconWidth*IconSize. Recommended: 2
  if ~isfield(config,'Menu'),             config.Menu             = {}; end % custom menu
  if ~isfield(config,'ExitConfirm'),      config.ExitConfirm      = 1; end  % Ask for confirmation before exiting iView if 'yes' or 1
  if ~isfield(config,'Version'),          config.Version          = '1.0'; end
  if ~isfield(config,'IconStyle'),        config.IconStyle        = 'togglebutton'; end % style of Icons: checkbox togglebutton. Recommended: togglebutton
  if ~isfield(config,'OutputFormat'),     config.OutputFormat     = 'pdf'; end  % prefered export format
  if ~isfield(config,'Labels'),           config.Labels           = {'Default','Background'}; end % predefined labels for data sets
  if ~isfield(config,'DragDropNoTooltips'),config.DragDropNoTooltips='yes'; end % 'yes' unactivates Tooltips, but activates DnD
  if ~isfield(config,'UseSystemDialogs'), config.UseSystemDialogs='no'; end % File selector specification: no= use iLoad selector, else defaults to Matlab built-in

function config=iView_private_config_save(config)
  % check config file name to use
  filename=fullfile(prefdir, 'iView.ini'); % store preferences in PrefDir (Matlab)
  config.FileName = filename;
  str = [ '% iView configuration script file ' sprintf('\n') ...
          '%' sprintf('\n') ...
          '% Matlab ' version ' m-file ' filename ' saved on ' datestr(now) ' with iview(''save_config'')' sprintf('\n') ...
          '%' sprintf('\n') ...
          '% The configuration must be specified as a structure with fields:' sprintf('\n') ...
          '%  PaperType        = ''A4'';                        % Paper type for iView window printing' sprintf('\n') ...
          '%  PaperUnits       = ''centimeters'';               % Paper unit specification for size/position' sprintf('\n') ...
          '%  PaperOrientation = ''landscape'';                 % Paper orientation' sprintf('\n') ...
          '%  PaperPosition    = [5 5 20 15];                 % Printing area on paper [Position Size ]' sprintf('\n') ...
          '%  Color            = [1 1 1];                     % Color of iView window background.' sprintf('\n') ...
          '%  IconSize         = 64;                          % Icon size in pixels' sprintf('\n') ...
          '%  IconWidth        = 2;                           % Real icon width is IconWidth*IconSize' sprintf('\n') ...
          '%  Menu             = {};                          % Custom menu items as { { ''Menu/Item'',''Callback''}, ... }' sprintf('\n') ...
          '%  ExitConfirm      = 1;                           % Ask for confirmation before exiting iView if ''yes'' or 1' sprintf('\n') ...
          '%  IconStyle        = ''togglebutton'';              % Style of Icons: checkbox togglebutton' sprintf('\n') ...
          '%  OutputFormat     = ''pdf'';                       % Prefered export format' sprintf('\n') ...
          '%  Labels           = {''Default'',''Background''};    % Predefined labels for data sets' sprintf('\n') ...
          '%  DragDropNoTooltips=''yes'';                       % ''yes'' unactivates Tooltips, but activates Drag and Drop' sprintf('\n') ...
          '%  UseSystemDialogs = ''no'';                        % File selector specification: no= use iLoad selector, else defaults to Matlab built-in' sprintf('\n') ...
          '%' sprintf('\n') ...
          '% NOTE: The resulting structure must be named "config"' sprintf('\n') ...
          '%' sprintf('\n') ...
          class2str('config', config) ];
  [fid, message]=fopen(filename,'w+');
  if fid == -1
    warning(['Error opening file ' filename ' to save iView configuration.' ]);
  else
    fprintf(fid, '%s', str);
    fclose(fid);
  end
  
