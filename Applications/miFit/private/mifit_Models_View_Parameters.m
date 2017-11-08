function stop=mifit_Models_View_Parameters(varargin)
  % [internal] mifit_Models_View_Parameters: the window showing Model parameters
  %   get 1st selected Model from Data set or Models menu current choice
  %   Display a uitable with columns:
  % [ Parameters | ParameterValues | ParameterUncertainty | constraints.fixed | constraints.min | constraints.max ]
  
  persistent cache

  if isempty(cache) || ~isstruct(cache)
    cache.numPars = 0;
    cache.modelTag='';
    cache.dataTag ='';
    cache.fig     =[];
    cache.table   =[];
  end
  
  stop=false;
  % 'fast' mode for a specific action, or update of the Parameters column.
  if nargin == 1 && ischar(varargin{1})
    action = varargin{1};
    stop = mifit_fig('mifit_View_Parameters');
    switch lower(action)
    case 'update'
      cache.modelTag = '';
    case 'histograms'
      mifit_Models_View_Parameters_Histograms;
      return
    end
  elseif nargin == 1 && numel(varargin{1}) == 1 && ishandle(varargin{1})
    stop = mifit_fig('mifit_View_Parameters');
    if ~isempty(stop), figure(stop); end
  elseif nargin == 1 && isnumeric(varargin{1}) && ~all(ishandle(varargin{1}))
    % quick update of the 'Value' column with the given parameter vector
    Parameters = varargin{1};
    f = mifit_fig('mifit_View_Parameters');
    if isempty(f) || ~ishandle(f), return; end
    t = getappdata(f, 'TableHandle');
    if ~ishandle(t), return; end
    Data = get(t, 'Data'); n = size(Data,1);
    if n == numel(Parameters)
      Data(:,2) = mat2cell(Parameters(:), ones(n,1),1);
      set(t, 'Data',          Data);
    end
    return
  end

  % get the Currently selected Data set
  if isappdata(mifit_fig, 'CurrentDataSet') && ~isempty(getappdata(mifit_fig, 'CurrentDataSet'))
    d              = getappdata(mifit_fig, 'CurrentDataSet');
    index_selected = getappdata(mifit_fig, 'CurrentDataSetIndex');
  else
    [d, index_selected] = mifit_List_Data_pull(); % get selected objects
    if numel(d)>1 && ~isempty(index_selected)
      d=d{1}; index_selected=index_selected(1);
    end
  end
  if iscell(d), d=d{1}; end
  
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

    t = uitable('Parent',f, ...
      'ColumnName',    options.ColumnName, ...
      'ColumnEditable',options.ColumnEditable, ...
      'ColumnFormat',  options.ColumnFormat, ...
      'FontSize',      config.FontSize, ...
      'RowStriping','on',...
      'CellEditCallback', @mifit_Models_View_Parameters_Edit, ...
      'Units','normalized', 'Position', [0 0 1 1]);

    % and trigger other updates (create all)
    resize                = 2;
    created               = true;
    data_or_model_changed = 1;
    % add contextual menu
    uicm = uicontextmenu('Parent',f); 
    uimenu(uicm, 'Label', [ 'Model:' model.Name ]);
    uimenu(uicm, 'Label', [ 'Data set:' char(d) ]);
    uimenu(uicm, 'Separator','on','Label','Copy all to clipboard', ...
      'Callback', @mifit_Models_View_Parameters_copy);
    uimenu(uicm, 'Label','Paste Values from clipboard', ...
      'Callback', @mifit_Models_View_Parameters_paste);
    uimenu(uicm, 'Label','Export to CSV file...', ...
      'Callback', @mifit_Models_View_Parameters_export);
    uimenu(uicm, 'Label','Plot parameter distributions...', ...
      'Callback', @mifit_Models_View_Parameters_Histograms);
    set(t,   'UIContextMenu', uicm);
    
    % we must store that window ID into the Appdata
    cache.table = t;
    setappdata(f, 'TableHandle',t);
  else
    t = getappdata(f, 'TableHandle');
    
    % different update levels:
    %   resize table when numel(Parameters) change
    resize = (n ~= cache.numPars);
    created               = false;
    %   update content when model.Tag or d.Tag has changed                 
    data_or_model_changed = (~isempty(d) && ~strcmp(cache.dataTag, d.Tag)) ...
      || ~strcmp(cache.modelTag, model.Tag);   
  end
  
  % the window must be resized as number of Parameters has changed
  if resize 
    % determine the new window size to show
    numCol    = 6;
    TextWidth = config.FontSize;
    TextHeight= config.FontSize*1.7;
    % height is given by the number of fields. Add 1 for the header line.
    height    = (n+1)*TextHeight;
    % width is given by the length of the longest RowName
    width     = ones(1,numCol+1)*10;
    width(1)  = max(cellfun(@numel,strtok(model.Parameters)));
    width     = width*TextWidth;
    options.ColumnWidth = mat2cell(width, 1, ones(1,numCol+1,1));
  
    % compare requested size to current window size
    p    = get(f, 'Position');
    if resize > 1, p(3) = sum(width); end
    p(4) = height;
    set(f, 'Position',p);
    % resize the window only at creation
    if created, set(t, 'ColumnWidth', options.ColumnWidth); end
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
        if ~isfield(modelValue.FitOutput.constraints, 'index_variable')
          modelValue.FitOutput.constraints.index_variable = 1:numel(dat.sig);
        end
        dat.sig(modelValue.FitOutput.constraints.index_variable) = ...
            max(dat.sig(modelValue.FitOutput.constraints.index_variable), ...
            modelValue.FitOutput.parsHessianUncertainty);
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
  
  stop = f;
  
% ------------------------------------------------------------------------------
function mifit_Models_View_Parameters_copy(varargin)
  t = getappdata(gcbf, 'TableHandle');
  d = get(t, 'Data');
  if isempty(d), return; end
  disp([ mfilename ': copy Parameters to the clipboard.' ]);
  c = [ get(gcbf, 'Name') sprintf('\n') ];
  for index=1:size(d,1)
    c = [ c sprintf('%20s %g %g\n', d{index,1:3}) ];
  end
  clipboard('copy', c);
  
function mifit_Models_View_Parameters_paste(varargin)
  str = clipboard('paste');
  t = getappdata(gcbf, 'TableHandle');
  d = get(t, 'Data');

  % check if this is a multiline string which can be evaluated as values (e.g. '0.01\n0.02 ...')
  val = str2num(str);
  if size(val,1) == size(d,1) % same number as lines
    for index=1:size(val,1); d{index,2} = val(index,1); end
    set(t, 'Data',          d);
  else
  % check if we find parameter names as first token
    str = textscan(str,'%s','delimiter','\n','whitespace',''); % read all lines
    str = str{1};
    % get par names from the Table
    tbl_names = strtok(d(:,1));
    for index=1:numel(str)
      if isempty(str{index}), continue; end
      str_name  = strtok(str{index});  % parameter name from the clipboard line
      str_val   = regexp(str{index},'\<(([1-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?\>','match');
      if isempty(str_val) || isempty(str_name), continue; end
      str_val=str2num(str_val{1});
      to_change = find(~cellfun(@isempty, strfind(tbl_names,str_name)));  % search it in the Table
      if isempty(to_change), continue; end  % not found in names
      d{to_change,2} = str_val(1);
    end
    set(t, 'Data',          d);
  end
  
function mifit_Models_View_Parameters_export(varargin)
  % exports the data set to a CSV file
  t = getappdata(gcbf, 'TableHandle');
  d = get(t, 'Data');
  if isempty(d), return; end
  filename = uiputfile(...
    {'*.csv','Comma Separated Value file (Excel, *.csv)';'*.*','All files'}, ...
    [ 'miFit:Parameters' ': Save as CSV (append to existing)' ]);
  if isempty(filename), return; end
  fid = fopen(filename,'a+');
  fprintf(fid, 'Model_Data,%s\n', get(gcbf,'Name'));
  fprintf(fid, 'Parameter,Value,Sigma_Half\n');
  for index=1:size(d,1)
    fprintf(fid, '%20s, %g, %g\n', d{index,1:3});
  end
  fclose(fid);
  
function mifit_Models_View_Parameters_Histograms(varargin)
  % display the Parameter distribution histograms
  f = mifit_fig('mifit_View_Parameters');
  if isempty(f) || ~ishandle(f), return; end
  o = getappdata(f, 'LastFitOutput');
  if isempty(o), return; end
  if iscell(o), o=o{1}; end
  figure('Name',[ 'Parameter distributions: ' datestr(now) ' ' o.optimizer ' ' o.model.Name ]);
  M=numel(o.parsBest);
  m=floor(sqrt(M)); n=ceil(M./m);
  index      = find(o.criteriaHistory < min(o.criteriaHistory)*10); % select a portion close to the optimum
  if numel(index)<20, index=1:numel(o.criteriaHistory); end
  
  % identify the parameters which are independent, from the Hessian Correlation matrix

  if ~isfield(o,'parsHessianCorrelation') || isempty(o.parsHessianCorrelation)
    corr = eye(M);
  else
    corr = o.parsHessianCorrelation;
  end
  parameters_independent = 1./sum(corr.^2);
  corr = corr - eye(M);

  % plot distributions
  for P=1:M
    subplot(m,n,P);
    [N,X]=hist(o.parsHistory(index,P));
    h=bar(X,N);
    if parameters_independent(P) < 0.7, 
      set(h, 'FaceColor','red');
      [~,most_correlated] = max(abs(corr(:,P)));
      if corr(most_correlated) < 0
        most_correlated = sprintf('\nCorrelated to -%s', o.parsNames{most_correlated});
      else
        most_correlated = sprintf('\nCorrelated to %s', o.parsNames{most_correlated});
      end
    else
      most_correlated = '';
    end
    title(sprintf('p(%i):%s = %g +/- %g%s', P, ...
      o.parsNames{P}, o.parsBest(P), ...
      o.parsHistoryUncertainty(P), most_correlated));
    xlabel(o.parsNames{P});
  end
