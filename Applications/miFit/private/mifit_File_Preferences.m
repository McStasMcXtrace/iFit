function config = mifit_File_Preferences(varargin)
% File/Preferences: open Preferences dialogue
% set directories to search for Models
% set FontSize (and update all Fonts in figure)
% set Save on exit
% save Preferences on dialogue close
  fig      = mifit_fig;
  config   = getappdata(mifit_fig, 'Preferences');
  % defaults are set when calling Preferences_Load at OpeningFcn (main/startup)
  filename = fullfile(prefdir, [ 'mifit' '.ini' ]);
  options.Name       = [ 'miFit' ': Preferences ' filename ];
  options.ListString = {'FontSize Font size [10-36]', ...
    'Models_Store Store Models when creation time is longer than [sec, 0:always, Inf:never, default=3]', ...
    'Models_Replot Replot Models after changing parameters when evauation time is shorter than [sec, 0:always, Inf:never, default=5]', ...
    'History_Level Undo levels to keep [2-50, reduce if you handle large/many data sets]', ...
    'Fit_Verbose Verbosity when Fitting [yes shows criteria, parameters and final distributions]', ...
    'ProxyHost Proxy address if you are behind a proxy [e.g. myproxy.mycompany.com or empty]', ...
    'ProxyPort Proxy port if you are behind a proxy [8888 or 0 or empty]' };
  options.FontSize   = config.FontSize;
  %options.Units      = 'normalized';
  %options.Position   = [ 0.02 0.02 0.98 0.98 ];
  options.ColumnEditable = true;
  options.TooltipString = sprintf([ 'Modify the %s Preferences.\n' ...
      'You can specify additional menu items by entering a cell (pairs)\n' ...
      '* Menu_<Label> = {''Item_Label'',''Command'', ...}' ], 'miFit');

  config1 = uitable(config, options); % structure GUI as a Table (spreadsheet)
  if isempty(config1), return; end
  config  = config1;
  % set new Preferences
  config.FontSize     =min(max(config.FontSize, 10),36);
  config.History_Level=min(max(config.History_Level, 2),50);
  % save the window size and location
  set(fig,'Units','pixels');
  config.Position = get(fig,'Position');
  
  setappdata(mifit_fig, 'Preferences', config);
  mifit_Preferences_Apply;
  mifit_Preferences_Save(config);
