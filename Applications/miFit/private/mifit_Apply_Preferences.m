function config=mifit_Apply_Preferences
  fig = mifit_fig;
  config = getappdata(fig, 'Preferences');
  
  % change all FontSize
  h=[ findobj(fig, 'Type','uicontrol') ; ...
      findobj(fig, 'Type','axes')    ; findobj(fig, 'Type','text') ; ...
      findobj(fig, 'Type','uipanel') ; findobj(fig, 'Type','uitable') ];
  set(h, 'FontSize', config.FontSize);
  
  % change fontsize in Parameter table
  f = mifit_fig('mifit_View_Parameters');
  if ~isempty(f)
    t = get(f, 'Children');
    set(t, 'FontSize', config.FontSize);
  end
  
  % handle User menus
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
      separator = 'off';
      flag_added_new_menu = true;
    else separator = 'on';  % to separate from already existing items
    end
    
    % add menu items as pairs {'Label','Command'}
    items = config.(field);
    for it=1:2:numel(items)
      item = items{it};
      cmd  = items{it+1};
      % check the menu item Label
      if ~ischar(item) || isempty(item)
        disp([ mfilename ': invalid Menu ' field ':item ' num2str(it) ' should be a char not ' class(item) '. Skipping.' ]);
      elseif ~ischar(cmd) && ~isa(cmd, 'function_handle')
        disp([ mfilename ': invalid Menu ' field ':' item ' callback. Should be a char of function_handle. Skipping.' ]);
      elseif ~isempty(item) && ~isempty(cmd)
        % check from separator as Label first char '|'
        if item(1) == '|', separator='on'; item(1)=[]; end
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
  % move the help menu to the right side
  menu_help  = findobj(fig, 'Tag','Menu_Help');
  right_help = copyobj(menu_help, fig);
  delete(menu_help);
  
  % for uimenu, could we use tip given by Y Altman 
  % <https://fr.mathworks.com/matlabcentral/newsreader/view_thread/148095>
