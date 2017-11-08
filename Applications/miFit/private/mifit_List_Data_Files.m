function mifit_List_Data_Files(varargin)
% [internal] mifit_List_Data_Files: called when clicking on the listbox
% TODO: support KeyPressFcn to e.g. 'p'=plot(log10(d)), 'f'=fit, 's'=save, 'o'=open, ...
% TODO: install a UIContextMenu in the Data sets list uicontrol ?
  if nargin && ishandle(varargin{1}), 
    obj=varargin{1}; 
  else obj=[]; end

  [d, index_selected] = mifit_List_Data_pull(); % get selected objects
  
  if isempty(index_selected) || numel(d) == 0
    setappdata(mifit_fig, 'CurrentDataSet', []);
    setappdata(mifit_fig, 'CurrentDataSetIndex', []);
  end

  if ~isempty(obj)
    if strcmp(get(gcbf,'SelectionType'), 'open')  % double-click
      mifit_Data_Plot([ d{:} ]);  % plot every time we change the data set selection
    end
  end
  
  % when a CurrentModel exists (previously selected from Models menu)
  % we assign that one to Data sets which have no Model defined
  if ~isempty(getappdata(mifit_fig, 'CurrentModel'))
    model = getappdata(mifit_fig, 'CurrentModel');
    D     = getappdata(mifit_fig, 'Data');  % all data sets
    model_assignments = 0;
    for index=index_selected(:)'
      this_d = D{index};
      if isempty(this_d) || ~isa(this_d,'iData'), continue; end
      previous_model = [];
      % get the Model stored in the Dataset (after fit)
      if isfield(this_d, 'Model')
        previous_model = get(this_d, 'Model');
      elseif ~isempty(findfield(this_d, 'Model'))
        previous_model = get(this_d, findfield(this_d, 'Model', 'cache first'));
      end
      if isempty(previous_model)  % store the CurrentModel in the Dataset
        this_d = setalias(this_d, 'Model', model);
        model_assignments = model_assignments+1;
      end
      D{index} = this_d;
    end
    if model_assignments
      setappdata(mifit_fig, 'Data', D);
      mifit_disp([ 'Assigning Model "' model.Name '" to ' num2str(model_assignments) ' Data set(s).' ]);
    end
  end
  
  % store the first selected Model so that its Model Parameters can be updated
  % in mifit_Models_View_Parameters
  if ~isempty(index_selected) && numel(d)
    d=d(1); index_selected=index_selected(1);
  end
  setappdata(mifit_fig, 'CurrentDataSet', d);
  setappdata(mifit_fig, 'CurrentDataSetIndex', index_selected);
  
  if ~isempty(mifit_fig('mifit_View_Parameters'))
    % trigger an update of the Parameter window when already opened
    mifit_Models_View_Parameters('update');
  end

