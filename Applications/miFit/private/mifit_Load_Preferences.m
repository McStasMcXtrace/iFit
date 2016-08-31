function config = mifit_Load_Preferences
  file = fullfile(prefdir, [ 'mifit' '.ini' ]);
  content = ''; config = '';
  if ~isempty(dir(file))
    try
      content = fileread(file);
      evalc(content);% this should make a 'config' variable
      mifit_disp([ '[Load_Preferences] Loading Preferences from ' file ]);
    end
  end
  if isempty(config)
    % default configuration
    config.FontSize         = max(12, get(0,'defaultUicontrolFontSize'));
    config.Save_Data_On_Exit= 'yes';
    config.Store_Models     = 3;  % time required for creation. Store when > 0:always, Inf=never
    config.History_Level    = 10;
  end
  setappdata(mifit_fig, 'Preferences', config);
  set(0,'defaultUicontrolFontSize', config.FontSize);
