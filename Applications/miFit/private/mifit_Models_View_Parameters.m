function mifit_Models_View_Parameters(varargin)
  % get 1st selected Model from Data set or Models menu current choice
  % Display a uitable with columns:
  % [ Parameters | ParameterValues | ParameterUncertainty | constraints.fixed | constraints.min | constraints.max ]
  
  persistent cache
  
  if isempty(cache) || ~isstruct(cache)
    cache.numPars = 0;
    cache.modelTag='';
    cache.dataTag ='';
    cache.fig     =[];
  end
  
  % get the Currently selected Data set
  if isappdata(mifit_fig, 'CurrentDataSet')
    d              = getappdata(mifit_fig, 'CurrentDataSet');
    index_selected = getappdata(mifit_fig, 'CurrentDataSetIndex');
  else
    [d, index_selected] = mifit_List_Data_pull(); % get selected objects
    if ~isempty(d) && ~isempty(index_selected)
      d=d(1); index_selected=index_selected(1);
    end
  end
  % get the Model to use.
  % if a Data set selection exists, get the first data set Model, or left empty
  model = []; modelValue = [];
  % search for a Model in the Data set
  if ~isempty(d)
    % get the Model stored in the Dataset (after fit)
    if isfield(d, 'Model')
      model = get(d, 'Model');
    elseif ~isempty(findfield(d, 'Model'))
      model = get(d, findfield(d, 'Model', 'cache first'));
    end
    % get the fit 'output' stored in the Dataset (after fit)
    if isfield(d, 'ModelValue')
      modelValue = get(d, 'ModelValue');
    elseif ~isempty(findfield(d, 'ModelValue'))
      modelValue = get(d, findfield(d, 'ModelValue', 'cache first'));
    end
  end
  % if no model from Data but CurrentModel exists, get it.
  if isempty(model)
    model = getappdata(mifit_fig, 'CurrentModel');
  elseif isa(model, 'iFunc')
    setappdata(mifit_fig, 'CurrentModel', model);
  end
  % else return (the User has to select a Data set or Model to view parameters)
  if isempty(model) || ~isa(model, 'iFunc')
    setappdata(mifit_fig, 'CurrentModel', []);
    return
  end
  
  % now display the figure and uitable
  % [ Parameters | ParameterValues | ParameterUncertainty | constraints.fixed | constraints.min | constraints.max ]   
  
  % window creation ------------------------------------------------------------
  % determine if the figure has to be created, or just updated

  % build the figure
  f = mifit_fig('mifit_View_Parameters');
  if isempty(f)
    if ishandle(cache.fig), f=cache.fig; end
  end
  n            = numel(model.Parameters);               % number of Parameters
  options.Name = [ 'mifit: ' model.Name ' ' char(d) ];  % Title
  config       = getappdata(mifit_fig, 'Preferences');  % config for FontSize
  
  % the window does not exist, we create it...
  if isempty(f)
    options.ColumnName    = {'Parameter','Value',  'Uncertainty','Fixed',  'Min',    'Max'    };
    options.ColumnFormat  = {'char',     'numeric','numeric',    'logical','numeric','numeric'};
    options.ColumnEditable= [ false       true      false         true      true      true    ];
    
    f = figure('Name',options.Name, 'MenuBar','none', ...
      'Tag','mifit_View_Parameters', ...
      'HandleVisibility','callback','NextPlot','new');

    % TODO: should install a callback when Parameters are changed: update the Data 
    % set Model Parameters, Constraints and replot.
    
    t = uitable('Parent',f, ...
      'ColumnName',    options.ColumnName, ...
      'ColumnEditable',options.ColumnEditable, ...
      'ColumnFormat',  options.ColumnFormat, ...
      'FontSize',      config.FontSize, ...
      'RowStriping','on',...
      'CellEditCallback', @mifit_Models_View_Parameters_Edit, ...
      'Units','normalized', 'Position', [0 0 1 1]);
    % and trigger other updates (create all)
    resize                = 1;
    data_or_model_changed = 1;
    % we must store that window ID into the Appdata
  else
    t = get(f, 'Children');
    
    % different update levels:
    %   resize table when numel(Parameters) change
    resize = (n ~= cache.numPars);
    %   update content when model.Tag or d.Tag has changed                 
    data_or_model_changed = (~isempty(d) && ~strcmp(cache.dataTag, d.Tag)) ...
      || ~strcmp(cache.modelTag, model.Tag);   
  end
  
  % the window must be resized as number of Parameters has changed
  if resize
    % determine the new window size to show
    numCol    = 6;
    TextWidth = config.FontSize*.8;
    TextHeight= config.FontSize*1.7;
    % height is given by the number of fields. Add 1 for the header line.
    height    = (n+1)*TextHeight;
    % width is given by the length of the longest RowName
    width     = ones(1,numCol+1)*5;
    width(1)  = max(cellfun(@numel,model.Parameters));
    width     = width*TextWidth;
    options.ColumnWidth = mat2cell(width, 1, ones(1,numCol+1,1));
  
    % compare requested size to current window size
    p    = get(f, 'Position');
    p(3) = sum(width);
    p(4) = height;
    set(f, 'Position',p);
    set(t, 'ColumnWidth', options.ColumnWidth);
  end
  
  % things to update -------------------------------------------------------------
  % determine if the Model/Data set is the same, in which case we only update the table Data
  if data_or_model_changed
    options.TooltipString = { model.Name, model.Description, char(d), ...
    'The Model will be updated when changing values.' };
    options.TooltipString = textwrap(options.TooltipString,80);
    options.TooltipString = sprintf('%s\n', options.TooltipString{:});
    set(t, 'TooltipString', options.TooltipString);
    set(f,'Name', options.Name);
  end
  
  % assemble the Data

  % after the fit, the iData object contains:
  % alias: Model
  % alias: ModelValue = output.modelValue with members:
  %   Constraints, FitOptions, FitOutput (history and uncertainties)
  % alias: ModelParameters
  
  if data_or_model_changed
    % get the Data content
    dat.pars = model.Parameters;
    
    dat.vals = model.ParameterValues;
    if numel(dat.vals) == n
      dat.vals = mat2cell(dat.vals(:), ones(n,1),1);
    else dat.vals = cell(n, 1); end
    
    if ~isempty(modelValue) && isfield(modelValue,'FitOutput')
      dat.sig = modelValue.FitOutput.parsHistoryUncertainty;
      if ~isempty(modelValue.FitOutput.parsHessianUncertainty)
        dat.sig = max(dat.sig, modelValue.FitOutput.parsHessianUncertainty);
      end
      dat.sig = mat2cell(dat.sig(:), ones(n,1),1);
    else
      dat.sig = cell(n, 1);
    end
    dat.min = mat2cell(model.constraint.min,   ones(n,1),1);
    dat.max = mat2cell(model.constraint.max,   ones(n,1),1);
    dat.fix = mat2cell(logical(model.constraint.fixed), ones(n,1),1);
    
    % change all NaN (might appear strange for dummy user) into empty
    for index=1:n
      if isnan(dat.min{index}),  dat.min{index}=[]; end
      if isnan(dat.max{index}),  dat.max{index}=[]; end
      if isnan(dat.vals{index}), dat.vals{index}=[]; end
    end
    % transfer it to the cell for uitable
    Data = cell(n,6);
    Data(:,1)= dat.pars;  % Parameter names:          model.Parameters
    Data(:,2)= dat.vals;  % Parameter values:         model.ParameterValues
    Data(:,3)= dat.sig;   % Parameter uncertainties:  output.parsHistoryUncertainty ...
    Data(:,4)= dat.fix;   % Parameter fixed/free:     model.constraint.fixed
    Data(:,5)= dat.min;   % Parameter min:            model.constraint.min
    Data(:,6)= dat.max;   % Parameter max:            model.constraint.max
    set(t, 'Data',          Data);
  end
 
% we store the Parameter window Tag so that consecutive calls can determine
% what has to be updated in consecutive calls.
  if ~isempty(d) cache.dataTag  = d.Tag; 
  else           cache.dataTag = ''; end
  cache.modelTag = model.Tag;
  cache.numPars  = n;
  cache.fig      = f;

