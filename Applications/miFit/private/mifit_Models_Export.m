function mifit_Models_Export(varargin)
% Models/Export: save user Models
  % get the list of 'static' iFunc models (which have been created and stored in the Models menu)
  [ifuncs, labels,indices] = mifit_Models_GetList();
  if isempty(indices), return; end

  % pop-up the iFunc.save export dialogue
  if numel(ifuncs) == 1
    save(ifuncs,'gui');
  else
    % pop-up a dialogue box to select those to export, with select all button
    [selection, ok] = listdlg('ListString', labels, 'SelectionMode', 'multiple', ...
      'Name','miFit: Select Models to Export', 'ListSize',[400 160]);
    if isempty(selection) || ok ~= 1, return; end
    save(ifuncs(selection),'gui');
  end
