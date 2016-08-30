function mifit_History_push
% append current Data into the History
  fig = mifit_fig;
  History = getappdata(fig, 'History');
  Data    = getappdata(fig, 'Data');
  if ~isempty(Data)
    History{end+1} = Data;
    if numel(History) > 10, History(1:(end-9)) = []; end
    setappdata(fig, 'History', History);
  end
