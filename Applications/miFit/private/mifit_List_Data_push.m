function mifit_List_Data_push(d, flag_replace)
% [internal] mifit_List_Data_push: put a new data set at the end of the stack
  if isempty(d),       return; end
  if nargin == 1, flag_replace = []; end
  if iscell(d)
    for index=1:numel(d)
      mifit_List_Data_push(d{index}, flag_replace);
    end
    return
  end
  if ~isa(d, 'iData')
    try
      mifit_disp(evalc('disp(d)'));
    end
    return; 
  end
  fig = mifit_fig;

  % update AppData Stack
  if numel(d) > 1, d = d(:); end
  Data = getappdata(fig, 'Data');

  if strcmp(flag_replace,'replace')
    % we search for data sets that have the same Tag, and replace them
    d0 = [];
    Tags = get(Data, 'Tag');
    for index=1:numel(d)
      match = find(strcmp(get(d(index),'Tag'), Tags));
      if ~isempty(match)
        if numel(Data) == 1, Data = d(index); 
        else Data(match) = d(index); % replace existing elements
        end
      else
        d0 = [ d0 ; d(index) ];      % append new elements
      end
    end
    d = d0;
  end
  
  Data = [ Data ; d ];  % a column of iData set
  setappdata(fig, 'Data', Data);
  
  % update the List labels by appending the Name at the end
  hObject        = mifit_fig('List_Data_Files');
  list           = get(hObject,'String');
  list0          = numel(list);
  index_selected = get(hObject,'Value');
  if max(index_selected) > numel(list), index_selected = []; end
  for index=1:numel(d)
      index_selected(end+1) = list0+index;
      list{end+1} = char(d(index));
  end
  set(hObject,'String', list, 'Value', index_selected);
  
  % Update the History with the new stack
  mifit_History_push;
  if numel(d) 
    mifit_disp('Importing into List:')
    mifit_disp(char(d))
  end
