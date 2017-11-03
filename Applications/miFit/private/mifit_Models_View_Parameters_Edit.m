function mifit_Models_View_Parameters_Edit(source, event)
% [internal] mifit_Models_View_Parameters_Edit:
%   callback when editing items in the Model Paremeter window

if nargin < 2, return; end
if isempty(source)   || isempty(event),   return; end
if ~ishandle(source) || (~isstruct(event) && ~isobject(event))
  return; end
if ~strcmpi(get(source,'Type'),'uitable'), return; end

% test if data changed, else return
if event.NewData == event.PreviousData, return; end

% we get the table content, and update the value
Data = get(source, 'Data');
Data{event.Indices(1),event.Indices(2)} = event.NewData;

% now store the new content into the Model. 

% search for the Model: is it in a uimenu or a Dataset ?
% That Model can be stored in the Data set and/or a uimenu item iFunc
% get that from the table.UserData

% if the uimenu is an expression, it must first be instantiated (create an iFunc).
% then we can not store the new change.

% the Handle to the DataSet is always: mifit_fig('List_Data_Files')
% when no selection, we should store the modified Model Parameters in 
%   getappdata(mifit_fig, 'CurrentModelHandle') when UserData holds an iFunc

d = []; model = [];
if ~isempty(getappdata(mifit_fig, 'CurrentDataSetIndex')) % Data set selection in List
  % get the single Model stored in the Dataset (after fit)
  d     = getappdata(mifit_fig, 'CurrentDataSet');
  if isfield(d, 'Model')
    model = get(d, 'Model');
  elseif ~isempty(findfield(d, 'Model'))
    model = get(d, findfield(d, 'Model', 'cache first'));
  end
elseif ~isempty(getappdata(mifit_fig, 'CurrentModel'))
  model = getappdata(mifit_fig, 'CurrentModel');
  handle= getappdata(mifit_fig, 'CurrentModelHandle');
end

if isempty(model), return; end

% check if the initial table has empty slots (which should then be converted to NaN or 0)
val = Data(:,event.Indices(2)); % a cell with the column where change took place
switch event.Indices(2)
case 4 % Fixed must be logical
  val(cellfun(@isempty, val)) = { false };
otherwise
  val(cellfun(@isempty, val)) = { NaN };
end
val = [ val{:} ];

% Event/update Signification depends on the Column index = event.Indices
switch event.Indices(2)
case 2
  % 2: Value -> store in Model.ParameterValues
  model.ParameterValues = val;
case 4
  % 4: Fixed -> store in Model.constraint.fixed
  model.Constraint.fixed = val;
case 5
  % 5: Min   -> store in Model.constraint.min
  model.Constraint.min = val;
case 6
  % 6: Max   -> store in Model.constraint.max
  model.Constraint.max = val;
otherwise
  disp([ mfilename ': ERROR: unsupported modified column ' num2str(event.Indices(2)) ]);
  return
end

if ~isempty(d)
  d = setalias(d, 'Model', model);
  % update the modelValue when the evaluation time is small (< .5 sec)
  if model.Duration > 0 && model.Duration < 0.5
    try
      % evaluate model with its parameters (Edit) and Data set axes
      if ndims(d) == ndims(model)
          modelValue = d(model);
          d = set(d, 'ModelValue', modelValue);
          % update plot, if found
          h = mifit_fig([ 'plot_' d.Tag ]);
          if ~isempty(h)
            plot(d,'light transparent grid tight replace');
          end
          % and put back focus to the Parameter Window
          figure(mifit_fig('mifit_View_Parameters'));
      end
    catch ME
      disp(getReport(ME))
    end
  end
  
  D     = getappdata(mifit_fig, 'Data');  % all data sets
  if numel(D) > 1 && getappdata(mifit_fig, 'CurrentDataSetIndex') <= numel(D)
    D{getappdata(mifit_fig, 'CurrentDataSetIndex')} = d;
    
    val = val(event.Indices(1));  % only the modified value is forwarded
    pars=model.Parameters;
    name= pars{event.Indices(1)};
    
    % now also set other selected data sets, when this is the same model
    % parameter values, fixed, limits
    index_selected = get(mifit_fig('List_Data_Files'),'Value');
    for index=index_selected(:)'
      % not the current data set (already updated)
      if getappdata(mifit_fig, 'CurrentDataSetIndex') == index, continue; end
      this_d = D{index};
      this_ = [];
      if isfield(this_d, 'Model')
        this_model = get(this_d, 'Model');
      elseif ~isempty(findfield(this_d, 'Model'))
        this_model = get(this_d, findfield(this_d, 'Model', 'cache first'));
      end
      if isempty(this_model), continue; end
      % the model must have similar parameters
      try
        switch event.Indices(2)
        case 2
          % 2: Value -> store in Model.ParameterValues
          this_model.(name) = val;
          this_model.(name) = val;  % duplicated in case no parameter values
                                    % were assigned before (none -> 0 -> val)
        case 4
          % 4: Fixed -> store in Model.constraint.fixed
          if isfinite(val) && val
            mlock(this_model,name);
          else
            munlock(this_model,name);
          end
        case 5
          % 5: Min   -> store in Model.constraint.min
          this_model.Constraint.min(event.Indices(1)) = val;
        case 6
          % 6: Max   -> store in Model.constraint.max
          this_model.Constraint.max(event.Indices(1)) = val;
        end
      end
      set(this_d, 'Model', this_model);
      D{index} = this_d;
    end  
  else
    D = { d };
  end
  setappdata(mifit_fig, 'Data', D);
  setappdata(mifit_fig, 'CurrentDataSet',d);
elseif ~isempty(handle), 
  set(handle, 'UserData',model);
end
