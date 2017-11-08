function mifit_File_Reset_Factory(ButtonName)
  % File Reset or factory settings
  if nargin == 0, ButtonName='Factory'; end
  if ~strcmpi(ButtonName, 'Reset') % Factory settings
    file = fullfile(prefdir, [ mfilename '.ini' ]);
    if ~isempty(dir(file)), delete(file); end
    file = fullfile(prefdir, [ mfilename '.mat' ]);
    if ~isempty(dir(file)), delete(file); end
    mifit_Preferences_Load();
    % mifit_Preferences_Save();
    if ~isempty(mifit_fig)
      setappdata(mifit_fig, 'Models',{});
      setappdata(mifit_fig, 'Optimizers',[]);
    end
  end
  if ~isempty(mifit_fig)
    hObject        = mifit_fig('List_Data_Files');
    set(hObject, 'String',[],'Value',[]);
    setappdata(mifit_fig, 'Data',{});
    setappdata(mifit_fig, 'CurrentDataSet', []);
    setappdata(mifit_fig, 'CurrentDataSetIndex', []);
    setappdata(mifit_fig, 'History', []);
    mifit_History_push();
  end

  file = fullfile(prefdir, [ mfilename '.log' ]);
  if ~isempty(dir(file)), delete(file); end
  
  if ~strcmpi(ButtonName, 'Reset') % Factory settings
    mifit_File_Exit;
    mifit;
  end
