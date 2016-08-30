function mifit_List_Data_push(d)
% put a new data set at the end of the stack
  if isempty(d),       return; end
  if ~isa(d, 'iData'), return; end
  fig = mifit_fig;

  % update AppData Stack
  if numel(d) > 1, d = d(:); end
  Data = getappdata(fig, 'Data');
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
  mifit_disp('Importing into List:')
  mifit_disp(char(d))
