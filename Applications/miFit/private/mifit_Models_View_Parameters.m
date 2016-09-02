function mifit_Models_View_Parameters(varargin)
  % get 1st selected Model from Data set or Models menu current choice
  % Display a uitable with columns:
  % [ Parameters | ParameterValues | ParameterUncertainty | constraints.fixed | constraints.min | constraints.max ]
  
  % get the Currently selected Data set
  if isappdata(mifit_fig, 'CurrentDataSet')
    d              = setappdata(mifit_fig, 'CurrentDataSet', d);
    index_selected = setappdata(mifit_fig, 'CurrentDataSetIndex', index_selected);
  else
    [d, index_selected] = mifit_List_Data_pull(); % get selected objects
    if ~isempty(d)
      d=d(1); index_selected=index_selected(1);
    end
  end
  % get the Model to use.
  % if a Data set selection exists, get the first data set Model, or left empty
  model = [];
  % search for a Model in the Data set
  if ~isempty(d)
    if isfield(d, 'Model')
      model = get(d, 'Model');
    elseif ~isempty(findfield(d, 'Model'))
      model = get(d, findfield(d, 'Model', 'cache first'));
    end
  end
  % if CurrentModel exists, get it.
  if isempty(model)
    model = getappdata(mifit_fig, 'CurrentModel');
  end
  % else return (the User has to select a Data set or Model to view parameters
  if isempty(model)
    return
  end
  
  % now display the figure and uitable
  % [ Parameters | ParameterValues | ParameterUncertainty | constraints.fixed | constraints.min | constraints.max ]
  
  % determine parameters for display -------------------------------------------
  
  
  n = numel(model.Parameters);
  
  
  % different update levels
  resize = 1;                     % when numel(Parameters) change
  data_or_model_changed = 1;
  model_changed = 1;
  
  % window creation ------------------------------------------------------------
  % determine if the figure has to be created, or just updated
  % if strcmp(d.Tag, '
  % build the figure
  f = mifit_fig('mifit_View_Parameters');
  
  % the window does not exist, we create it...
  if isempty(f)
    options.ColumnName = {'Parameter','Value','Uncertainty','Fixed','Min','Max'};
    options.ColumnFormat={'char',     'numeric','numeric',  'logical','numeric','numeric'};
    options.ColumnEditable=[false true false true true true ];
    config = getappdata(mifit_fig, 'Preferences');
    
    f = figure('Name',options.Name, 'MenuBar','none', 'Tag','mifit_View_Parameters');

    % TODO: should install a callback when Parameters are changed: update the Data 
    % set Model Parameters and replot.
    
    t = uitable('Parent',f, ...
      'ColumnName',    options.ColumnName, ...
      'ColumnEditable',options.ColumnEditable, ...
      'ColumnFormat',  options.ColumnFormat, ...
      'FontSize',      config.FontSize, ...
      'RowStriping','on',...
      'Units','normalized', 'Position', [0 0 1 1]);
  elseif  
    t = get(f, 'Children');
  end

  % determine if the Model/Data set is the same, in which case we only update the table Data
  
  if resize
    % determine the window size to show
    numCol = 6;
    TextWidth = options.FontSize*.8;
    TextHeight= options.FontSize*1.7;
    % height is given by the number of fields. Add 1 for the header line.
    height = (n+1)*TextHeight;
    % width is given by the length of the longest RowName
    width    = ones(1,numCol+1)*5;
    width(1) = max(cellfun(@numel,model.Parameters));
    width = width*TextWidth;
    options.ColumnWidth = mat2cell(width, 1, ones(1,numCol+1,1));
  
    % compare requested size to current window size
    p = get(f, 'Position');
    p(3) = sum(width);
    if p(4) > height, p(4) = height; end
    set(f, 'Position',p);
    set(t, 'ColumnWidth',   options.ColumnWidth);
  end
  
  % things to update -------------------------------------------------------------
  

  if data_or_model_changed
    options.TooltipString = { model.Description, char(d), ...
    'The Model will be updated when changing values.' };
    options.TooltipString = textwrap(options.TooltipString,80);
    options.TooltipString =sprintf('%s\n', options.TooltipString{:});
    set(t, 'TooltipString', options.TooltipString);
  end
  if model_changed
    options.Name       = [ 'mifit: ' model.Name ];
    set(f,'Name', options.Name);
  end
  
  % assemble the Data
  e = cell(n, 1);   % empty
  TF= mat2cell(false(n,1), ones(n,1),1);
  z = mat2cell(zeros(n,1), ones(n,1),1);
  mp = model.Parameters;
  mpv= model.ParameterValues;
  if isempty(mpv), mpv=z; 
  else mpv = mat2cell(mpv, ones(n,1),1); end
  Data = cell(n,6);
  Data(:,1)= mp;    % Parameter names
  Data(:,2)= mpv;   % Parameter values
  Data(:,3)= z;     % Parameter uncertainties
  Data(:,4)= TF;    % Parameter fixed/free
  Data(:,5)= e;     % Parameter min
  Data(:,6)= e;     % Parameter max
  
  set(t, ...
    'Data',          Data);
 

% we store the Parameter window Tag so that consecutive calls can determine
% what has to be updated in consecutive calls.

