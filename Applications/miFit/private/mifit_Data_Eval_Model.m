function mifit_Data_Eval_Model(varargin)
% Data/Evaluate Model
  if nargin && isa(varargin{1}, 'iData'), d=varargin{1}; 
  else d = mifit_List_Data_pull; end
  if all(isempty(d)), return; end
  set(mifit_fig,'Pointer','watch');
  model = {};
  if numel(d) == 1
    if ~iscell(d), d={ d }; end
    if isfield(d{1}, 'Model')
      model = { get(d{1}, 'Model') };
    elseif ~isempty(findfield(d{1}, 'Model'))
      model = { get([ d{1} ], findfield(d{1}, 'Model', 'cache first')) };
    end
  else
    model = getalias([ d{:} ], 'Model');
  end
  mifit_disp([ '[Data Eval Model] evaluating Model value for data sets ...' ]);
  for index=1:numel(d)
    this = d{index};
    if iscell(model) && index <= numel(model) && ~isempty(model{index})
      this_model = model{index};
      modelValue = this(this_model);  % eval iFunc
      this = set(this, 'ModelValue', modelValue);
      this = set(this, 'Model',      modelValue.Model);
      setappdata(mifit_fig, 'CurrentModel', modelValue.Model);
      mifit_disp(char(this))
      mifit_disp([ '  [ ' num2str(modelValue.Model.ParameterValues(:)') ' ]' ]);
      d{index} = this;
    end
  end
  mifit_List_Data_push(d, 'replace');  % replace existing data sets
  set(mifit_fig,'Pointer','arrow');
  
  if ~isempty(mifit_fig('mifit_View_Parameters'))
    % trigger an update of the Parameter window when already opened
    mifit_Models_View_Parameters('update');
  end
  
  mifit_History_push;
