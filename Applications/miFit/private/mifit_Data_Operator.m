function mifit_Data_Operator(op)
% Data/Transform: apply a binary operator on Datasets

  % we display a dialogue to select which operator to apply on the selected Datasets
  [d, index_selected] = mifit_List_Data_pull;
  if isempty(index_selected) || numel(d) == 0, return; end
  
  op = sort(op); % alpha order sort
    
  % show operator selection dialogue
  PromptString = cellstr(char([ d{:} ]));
  if numel(PromptString) > 3
    PromptString=[ PromptString(1:3) ; ' ...' ];
  end
  PromptString = [ PromptString ; 'Refer to the Help:Math Operators for more.' ];
  [selection, ok] = listdlg('ListString', op, 'ListSize',[500 160], ...
    'Name',[ 'miFit: Select operator to apply on ' ...
      num2str(numel(index_selected)) ' data set(s)' ], ...
      'PromptString', PromptString, 'SelectionMode','single');
  if isempty(selection) || ok ~= 1, return; end
  op = op{selection};
 
  % apply operator on Datasets
  mifit_disp([ 'Applying operator "' op '" on data set(s):' ]);
  mifit_disp(char([ d{:} ]))
  set(mifit_fig,'Pointer','watch');
  
  try
    [op, arg] = strtok(op);
    arg = str2num(strtok(arg));
    if isempty(arg)
      d = feval(op, [ d{:} ]);
    else
      % special case for operators with 2nd numeric argument
      d = feval(op, arg, [ d{:} ]);
    end
  catch ME
    mifit_disp([ '[Data operation] Failed applying operator "' op '"' ]);
    mifit_disp(getReport(ME));
    return
  end
    
  % upload new Data sets
  mifit_List_Data_push(d);
  set(mifit_fig,'Pointer','arrow');
