function mifit_Edit_Select_All(hObject, select)
% Edit/Select All/None: set the List selected values to all ones
% the select argument can be false or true to deselect/select all
  hObject = mifit_fig('List_Data_Files');
  items   = get(hObject,'String');
  index_selected = get(hObject,'Value');
  if nargin < 2, select=[]; end
  if (numel(index_selected) == numel(items) && isempty(select)) ...
    || (~isempty(select) && islogical(select) && ~select), index_selected = [];
  elseif isempty(select) || (~isempty(select) && islogical(select) && select), index_selected=1:numel(items); 
  elseif isnumeric(select) && min(select) > 0 && max(select) <= numel(items)
    index_selected = select;
  end
  if ~isempty(index_selected) && (min(index_selected) <= 0 || max(index_selected) > numel(items))
    index_selected = [];
  end
  set(hObject,'Value', index_selected);
  D     = getappdata(mifit_fig, 'Data');
  if numel(D) && ~isempty(index_selected)
    setappdata(mifit_fig, 'CurrentDataSet',    D{index_selected(1)});
    setappdata(mifit_fig, 'CurrentDataSetIndex', index_selected(1));
  else
    setappdata(mifit_fig, 'CurrentDataSet', []);
    setappdata(mifit_fig, 'CurrentDataSetIndex', []);
  end
