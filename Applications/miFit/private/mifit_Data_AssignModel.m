function mifit_Data_AssignModel(varargin)
% [internal] mifit_Models_GetList: this function is called when the user selects a Model from the Models menu
%   as a Callback
  
  model = get(varargin{1},'UserData'); % an iFunc or char/cellstr stored into UserData of the menu item.
  handle= [];
  
  if iscellstr(model), model = char(model)'; end
  if ischar(model)
    config = getappdata(mifit_fig, 'Preferences');
    try
      % create the Model from the stored expression
      tstart   = tic;
      modelF   = eval(model);
      telapsed = toc(tstart);
      if ~isempty(modelF) 
        if telapsed > config.Store_Models
          % push Model into the Models menu (to the Deck/Library/Stack)
          % all models are saved when exiting (when set in Preferences)
          handle = mifit_Models_Add_Entry(modelF);
          mifit_disp([ '[Data_AssignModel] Adding Models/' num2str(modelF.Dimension) ' menu item ' model ' ' modelF.Name ' [' modelF.Tag ']' ]);
        end
        model = modelF;
      end
    catch
      mifit_disp([ '[Data_AssignModel] Invalid Model expression ' model '. Skipping.' ]);
      return
    end
  elseif isa(model, 'iFunc')
    handle = varargin{1};
    model  = copyobj(model); % to get a new ID
  end

  if isempty(model) || ~isa(model, 'iFunc')
    setappdata(mifit_fig, 'CurrentModel',       []);
    setappdata(mifit_fig, 'CurrentModelHandle', []);
    return
  end
  setappdata(mifit_fig, 'CurrentModel',       model);  % store current selected Model and its Handle
  setappdata(mifit_fig, 'CurrentModelHandle', handle);
  
  % get selected Data sets indices in List and store the Model in there
  index_selected = get(mifit_fig('List_Data_Files'),'Value');
  D = getappdata(mifit_fig, 'Data');  % all data sets
  if numel(D) == 0 || isempty(index_selected), 
    mifit_disp([ 'Selected Model "' model.Name '".' ]);
    figure; plot(model);
    return; 
  end
  mifit_disp([ 'Assigning Model "' model.Name '" to ' num2str(numel(index_selected)) ' Data set(s):' ]);
  mifit_History_push();
  for index=index_selected(:)'
    if numel(D) > 1, this_d = D(index); else this_d = D; end
    this_d = setalias(this_d, 'Model', model);
    mifit_disp(char(this_d));
    if numel(D) > 1, D(index) = this_d; else D = this_d; end
    if index==index_selected(1), 
      setappdata(mifit_fig, 'CurrentDataSet', this_d); 
    end
  end
  setappdata(mifit_fig, 'Data', D);
  
  % trigger an update of the Parameter window when already opened
  if ~isempty(mifit_fig('mifit_View_Parameters'))
    mifit_Models_View_Parameters();
  end
