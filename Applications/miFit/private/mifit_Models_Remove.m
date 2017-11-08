function mifit_Models_Remove(varargin)
% Models/Remove: clear user models
  % select 'static' iFunc models to remove from the menu
   
  % get the list of 'static' iFunc models (which have been created and stored in the Models menu)
  [ifuncs, labels, indices,handles] = mifit_Models_GetList();
  if isempty(indices), return; end
  % pop-up a dialogue box to select those to remove, with select all button
  [selection, ok] = listdlg('ListString', labels, 'SelectionMode', 'multiple', ...
    'Name','miFit: Select Models to Remove', 'ListSize',[400 160]);
  if isempty(selection) || ok ~= 1, return; end
  delete(handles(selection))
  models = getappdata(mifit_fig,'Models');
  models(indices(selection)) = [];
  setappdata(mifit_fig,'Models',models);
