function config = mifit_Load_Preferences(file)
  if nargin == 0, file = []; end
  if isempty(file)
    file = fullfile(prefdir, [ 'mifit' '.ini' ]);
  end
  content = ''; config = '';
  if ~isempty(dir(file))
    try
      content = fileread(file);
      evalc(content);% this should make a 'config' variable
      mifit_disp([ '[Load_Preferences] Loading Preferences from ' file ]);
    end
  end
  if isempty(config) && nargin == 0
    % default configuration
    config.FontSize         = max(12, get(0,'defaultUicontrolFontSize'));
    config.Save_Data_On_Exit= 'yes';
    config.Store_Models     = 3;  % time required for creation. Store when > 0:always, Inf=never
    config.History_Level    = 10;
  end
  if isstruct(config)
    setappdata(mifit_fig, 'Preferences', config);
    if isfield(config, 'FontSize')
      set(0,'defaultUicontrolFontSize', config.FontSize);
    end
  end
  
  % resize main panel
  handle = findobj(mifit_fig, 'Tag','Panel_DataSets');
  set(handle, 'Units','normalized', 'Position',[0.002 0 0.99 0.99]);
