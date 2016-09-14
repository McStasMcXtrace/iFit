function mifit_Data_AssignModel(varargin)
  % this function is called when the user selects a Model from the Models menu.
  
  model = get(varargin{1},'UserData'); % an iFunc or char/cellstr stored into UserData of the menu item.
  
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
          mifit_Models_Add_Entry(modelF);
          mifit_disp([ '[Data_AssignModel] Adding Models/' num2str(modelF.Dimension) ' menu item ' model ' ' modelF.Name ' [' modelF.Tag ']' ]);
        end
        model = modelF;
      end
    catch
      mifit_disp([ '[Data_AssignModel] Invalid Model expression ' model '. Skipping.' ]);
      return
    end
  elseif isa(model, 'iFunc')
    model = copyobj(model); % to get a new ID
  end
  if isempty(model), return; end
  setappdata(mifit_fig, 'CurrentModel', model);  % store current selected Model
  % get selected Data sets indices in List
  index_selected = get(mifit_fig('List_Data_Files'),'Value');
  D = getappdata(mifit_fig, 'Data');  % all data sets
  if numel(D) == 0 || isempty(index_selected), 
    mifit_disp([ 'Selected Model "' model.Name '".' ]);
    figure; plot(model);
    return; 
  end
  mifit_disp([ 'Assigning Model "' model.Name '" to ' num2str(numel(index_selected)) ' Data set(s):' ]);
  mifit_History_push();
  if numel(D) > 1
    for index=index_selected(:)'
      D(index) = setalias(D(index), 'Model', model);
      mifit_disp(char(D(index)));
    end
  else
    D = setalias(D, 'Model', model);
    mifit_disp(char(D));
  end
  setappdata(mifit_fig, 'Data', D);
