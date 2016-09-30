function [d, index_selected]=mifit_List_Data_pull(varargin)
% get the selected Data List
% return the selected objects
  hObject = mifit_fig('List_Data_Files');
  d = [];

  index_selected = get(hObject,'Value');
  if isempty(index_selected), return; end
  
  d   = getappdata(mifit_fig, 'Data');
  if numel(d) > 1
      d = d(index_selected);
  end
