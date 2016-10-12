function mifit_History_pull
% [internal] mifit_History_pull for Edit/Undo: restore the last History state and deletes it
  fig = mifit_fig;
  History = getappdata(fig, 'History');
  if isempty(History), return; end
  Data         = History{end};
  History(end) = [];
  setappdata(fig, 'History', History);
  setappdata(fig, 'Data',    Data);
