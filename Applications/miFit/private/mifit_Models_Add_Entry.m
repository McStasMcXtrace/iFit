function mifit_Models_Add_Entry(model)
  % add a new Model (iFunc or cellstr) into the Models sub-menus.
  % the input argument should be a structure with members:
  %   model.callback
  %   model.label
  %   model.object
  % or an iFunc array
  % or a cellstr (expressions to evaluate at callback)
  
  % usage:
  % mifit_Models_Add_iFunc(gauss)         - must add a Gaussian entry in sub-menu, with callback to 'copyobj(gauss)'
  % mifit_Models_Add_iFunc('gauss+lorz')  - must add a 'gauss+lorz' entry in sub-menu, with callback to 'gauss+lorz'
  
  if nargin == 0 || all(isempty(model)), return; end
  
  % check different type of input arguments Models
  if ischar(model), model = cellstr(model); end % to avoid using numel on char
  
  % handle array of entries
  if numel(model) > 1  % model is a cell or array of entries
    mifit_disp([ '[Models] Setting ' num2str(numel(model)) ' Models...' ]);
    for index=1:numel(model)
      if iscell(model), mifit_Models_Add_Entry(model{index});
      else              mifit_Models_Add_Entry(model(index));
      end
    end
    return
  end
  
  % handle single entry
  % determine
  %   label: a string to display for the menu entry
  %   callback: a char/cellstr (expression) or iFunc model, stored in the menu 
  %     item UserData for assignement when selected
  %   dim:   dimensionality (to identify the sub-menu)
  
  label = ''; dim = []; callback = [];
  if all(isempty(model)), return; end
  
  % determine the 'callback' for the menu item
  if iscellstr(model)
    callback = char(model)';
  elseif isstruct(model) && isfield(model,'callback')
    callback = model.callback;
    if isfield(model, 'label'),     label = model.label; end
    if isfield(model, 'object'),    model = model.object; 
    elseif isfield(model, 'model'), model = model.model; 
    end
  elseif ischar(model)
    callback = model;
  elseif ~isa(model, 'iFunc')
    mifit_disp([ '[Models_Add_Entry] Invalid Model type ' class(model) '. Must be iFunc/cell/char/struct with "callback" member. Skipping.' ]);
  end
  
  % determine the dimensionality
  if isa(model, 'iFunc') && ~isempty(model)
    dim      = model.Dimension;
    if isempty(callback) callback = model; end
  else model = [];
  end

  if isa(callback,'iFunc') && ~isempty(callback)
    if isempty(dim), dim = callback.Dimension; end
    if isempty(label)
      label    = [ callback.Name ' [' callback.Tag ']' ];
    end
  end
  
  if isempty(dim) % this is a model creator (new instance) from expression
    % we assume this is an expression and try to evaluate it as such
    try
      modelF    = feval(callback, 'identify');
    catch
      mifit_disp([ '[Models_Add_Entry] Invalid Model expression ' callback '. Skipping.' ]);
      return
    end
    if isempty(modelF), return; end
    dim   = modelF.Dimension;
    if isempty(label), label = [ '"' callback '" = ' modelF.Name ]; end
  end
  % determine the name of the sub-menu to use
  if isempty(dim), dim = -1; end  % will use Others
  if dim > 0, model_submenu_name = sprintf('%dD', dim);
  else        model_submenu_name = 'Others'; end
  
  % create the sub-menu, if missing
  hmodels = mifit_fig('Menu_Model');
  submenu_handle = findobj(hmodels, 'Type','uimenu', 'Label', model_submenu_name);
  if isempty(submenu_handle)
    submenu_handle = uimenu(hmodels, 'Label', model_submenu_name);
  end
  
  % check if the model entry already exists
  children = findobj(submenu_handle, 'Label', label);
  if isempty(children)
    children = findobj(submenu_handle, 'UserData', callback);
  end
  if ~isempty(children)
    mifit_disp([ '[Models_Add_Entry] ' label ' is already in the list of usable Models. Skipping.' ]);
    return; 
  end
    
  handle = uimenu(submenu_handle, 'Label', label, 'UserData', callback, ...
                'CallBack', 'mifit(''Data_AssignModel'',gcbo)');
                
  % store the entry in the appdata
  % Models is a cell of entries which can be:
  %   struct.callback
  %   struct.callback and struct.label
  %   iFunc
  %   char (expression)
  Models = getappdata(mifit_fig, 'Models');
  if ~isa(model, 'iFunc'), model=[]; end
  Models{end+1} = struct('callback',callback, 'label',label, ...
    'object',model, 'handle',handle, 'dimension',dim);
  setappdata(mifit_fig, 'Models',Models);
  
  
