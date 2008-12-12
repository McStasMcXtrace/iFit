function config=iView_private_config(instance, action)

switch action
case 'load'
  config=iView_private_config_load(instance);
otherwise
  disp([' Unknown action ' action ' in ' mfilename ]);
end

% =============================================================================
%                    INLINE functions
% =============================================================================
function config=iView_private_config_load(instance)
config=[];

% is there an existing config ? If yes, re-use it
if ishandle(instance)
  if isappdata(instance, 'Config')
    config=getappdata(instance, 'Config');
  end
else % iView instance does not exist, search for an other instance to get config
  % find all existing iView instances
  instance_list=findall(0,'Tag','iView_instance');
  if ~isempty(instance_list)
    instance=instance_list(1);  % use first existing instance
    config=getappdata(instance, 'Config');
  end
end

if isempty(config) % doen when building first iView instance
  % else get default default configuration
  config.FileName         = 'iView';
  config.PaperUnits       = 'centimeters';
  config.PaperOrientation = 'landscape';
  config.PaperPosition    = [5 5 20 15];
  config.Color            = [1 1 1]; % white
  config.IconSize         = 64;

  % customized menus 
  % Label=cell of {'Menu/Item','callback'...} with associated callbacks
  % To insert a Separator, use 'Menu/Item'='Menu/Separator'
  config.Menu= {
  'Tools/Custom entry','disp(''hello'')', ...
  };

  % and overrides by any local configuration file
  
end

