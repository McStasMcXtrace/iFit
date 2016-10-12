function mifit_History_push
% [internal] mifit_History_push: append current Data into the History
  fig = mifit_fig;
  History = getappdata(fig, 'History');
  Data    = getappdata(fig, 'Data');
  config  = getappdata(fig, 'Preferences');
  if ~isfield(config, 'History_Level'), config.History_Level=10; end
  if ~isempty(Data)
    History{end+1} = Data;
    if numel(History) > config.History_Level, History(1:(end-config.History_Level+1)) = []; end
    setappdata(fig, 'History', History);
  end
