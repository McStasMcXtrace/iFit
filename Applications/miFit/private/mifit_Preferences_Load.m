function config = mifit_Preferences_Load(file)
% File/Load Preferences
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
    catch
      mifit_disp([ '[Load_Preferences] Error loading ' file ]);
    end
  end
  if isempty(config) && nargin == 0
    % request default configuration
    config = struct();
  end
  if isstruct(config)
    % make sure we have defined all defaults
    if ~isfield(config, 'FontSize'),          config.FontSize         = max(12, get(0,'defaultUicontrolFontSize')); end 
    if ~isfield(config, 'Models_Store'),      config.Models_Store     = 3; end
    if ~isfield(config, 'Models_Replot'),     config.Models_Replot    = 5; end
    if ~isfield(config, 'History_Level'),     config.History_Level    = 10; end
    if ~isfield(config, 'Fit_Verbose'),       config.Fit_Verbose      = 'no'; end
    if ~isfield(config, 'ProxyHost'),         config.ProxyHost        = ''; end
    if ~isfield(config, 'ProxyPort'),         config.ProxyPort        = 0; end
    if ~isempty(mifit_fig), setappdata(mifit_fig, 'Preferences', config); end
    if isfield(config, 'FontSize')
      set(0,'defaultUicontrolFontSize', config.FontSize);
    end
  end
  
  % resize main panel.
  if 0 && ~isempty(mifit_fig)
    handle = findobj(mifit_fig, 'Tag','Panel_DataSets');
    set(handle, 'Units','normalized', 'Position',[0.002 0 0.99 0.99]);
  end
