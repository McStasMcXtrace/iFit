function config=mifit_Preferences_Apply(config)
% [internal] mifit_Preferences_Apply: apply Preferences after load, from default 
%   or specified config structure
%   This function also handles User menus from the Preferences

  % get the Preferences
  fig = mifit_fig;
  if nargin == 0, config = []; end
  if ~isstruct(config) 
    config0 = getappdata(fig, 'Preferences'); 
    % config=struct: merge config and config0
    if isstruct(config)
      for f=fieldnames(config0)'
        if ~isfield(config, f{1}) config.(f{1}) = config0.(f{1}); end
      end
    else
      config = config0;
    end
  end
  
  % change all FontSize
  h=[ findobj(fig, 'Type','uicontrol') ; ...
      findobj(fig, 'Type','axes')    ; findobj(fig, 'Type','text') ; ...
      findobj(fig, 'Type','uipanel') ; findobj(fig, 'Type','uitable') ];
  if isfield(config, 'FontSize')
    set(h, 'FontSize', config.FontSize);
  end
  
  % change fontsize in Parameter table
  f = mifit_fig('mifit_View_Parameters');
  if ~isempty(f) && isfield(config, 'FontSize')
    t = findobj(f, 'Type','uitable');
    set(t, 'FontSize', config.FontSize);
  end
  
  % change window size/location
  Position = [];
  if isfield(config, 'Position')
    Position = config.Position;
  end
  if numel(Position) == 4
    set(fig, 'Units','pixels','Position',Position);
  elseif numel(Position) == 2 % size only
    set(fig, 'Units','pixels');
    p0 = get(fig,'Position');
    Position(1:2) = p0(1:2);
    set(fig, 'Position',Position);
  end
  
  % handle Proxy settings
  if isfield(config, 'ProxyHost') 
    if ~isempty(config.ProxyHost) && ischar(config.ProxyHost)
      java.lang.System.setProperty('http.proxyHost',config.ProxyHost);  % old Matlab stuff
      try
        com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxy(true);
        com.mathworks.mlwidgets.html.HTMLPrefs.setProxyHost(config.ProxyHost);
      end
      % proxy port
      if ~isempty(config.ProxyPort)
        if isnumeric(config.ProxyPort) && isscalar(config.ProxyPort) && config.ProxyPort > 0
          config.ProxyPort = num2str(config.ProxyPort); 
        elseif ~ischar(config.ProxyPort)
          config.ProxyPort = '';
        end
        java.lang.System.setProperty('http.proxyPort',config.ProxyPort);
        try
          com.mathworks.mlwidgets.html.HTMLPrefs.setProxyPort(config.ProxyPort);
        end
      end
    else  % no proxy
      java.lang.System.setProperty('http.proxyHost','');
      java.lang.System.setProperty('http.proxyPort','');
      try
        com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxy(false);
      end
    end
  end
    
  
  % handle User menus ==========================================================
  % dynamically build any other 'Menu_' items in the configuration
  % example: config.Menu_Tools={'Label','command'}; will add the label in the Menu
  flag_added_new_menu = false;
  for f=fieldnames(config)'
    % look for the Menu where to add items
    field = f{1}; menu_label='';
    if ~iscell(config.(field)),           continue; end % not a cell
    if mod(numel(config.(field)),2) ~= 0, continue; end % not even numel
    % get the desired target Menu
    if strncmp(lower(field), 'menu_', 5) % a new menu entry in config
      % get the Menu label
      menu_label = lower(field(6:end)); menu_label(1) = upper(menu_label(1));
      menu_label = strrep(menu_label, '_',' ');
    end
    if isempty(menu_label), continue; end
    
    % identify the menu where to add items (or create it)
    menu_handle = findobj(fig, 'tag', [ 'Menu_' genvarname(menu_label) ]);
    if isempty(menu_handle), menu_handle = findobj(fig, 'Label', menu_label); end
    % test if the Menu handle (when exists) is a root menu (no Callback)
    if ~isempty(get(menu_handle, 'CallBack')) menu_handle=[]; end
    % when the menu label does not exist yet as a Menu, we create it...
    if isempty(menu_handle)
      mifit_disp([ '[Init] Adding menu "' menu_label '" from Preferences...' ]);
      menu_handle = uimenu(fig, 'Label', menu_label, 'Tag', [ 'Menu_' genvarname(menu_label) ]);
      flag_added_new_menu = true;
    end
    
    % add menu items as pairs {'Label','Command'}
    items = config.(field);
    for it=1:2:numel(items)
      item = items{it};
      cmd  = items{it+1};
      % check the menu item Label
      if ~ischar(item) || isempty(item)
      elseif ~ischar(cmd) && ~isa(cmd, 'function_handle')
        disp([ mfilename ': invalid Menu ' field ':' item ' callback. Should be a char of function_handle. Skipping.' ]);
      elseif ~isempty(item) && ~isempty(cmd)
        % check from separator as Label first char '|'
        if item(1) == '|', separator='on'; item(1)=[]; else separator = 'off'; end
        % look if the menu item already exists
        item_handle = findobj(fig, 'Tag', [ menu_label '_' genvarname(item) ]);
        if isempty(item_handle)
          item_handle = findobj(fig, 'Label', item);
        end
        if ~isempty(item_handle)
          % already there: we update the Callback
          set(item_handle, 'Label', item, 'Callback', cmd, 'Tag', [ menu_label '_' genvarname(item) ], 'Separator', separator);
        else 
          % add the new menu entry
          item_handle = uimenu(menu_handle, ...
            'Label', item, 'Callback', cmd, 'Tag', [ menu_label '_' genvarname(item) ], 'Separator', separator);
        end
        separator = 'off';
      end
    end
  end
  
  % for uimenu, could we use tip given by Y Altman 
  % <https://fr.mathworks.com/matlabcentral/newsreader/view_thread/148095>
