function [d, index_selected]=mifit_List_Data_pull(varargin)
% [internal] mifit_List_Data_pull: get the selected Data List
% return the selected objects
  hObject = mifit_fig('List_Data_Files');
  d = [];

  index_selected = get(hObject,'Value');
  if isempty(index_selected), return; end
  
  index_selected = sort(index_selected);
  d   = getappdata(mifit_fig, 'Data');

  if numel(d) > 1
      % index_selected must only have items up to numel(d)
      index_selected = index_selected(index_selected <= numel(d));
      d = d(index_selected);
  end
