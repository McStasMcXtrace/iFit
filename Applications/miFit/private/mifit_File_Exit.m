function mifit_File_Exit(varargin)
% File/Exit: Quit and Save Data
  config = getappdata(mifit_fig, 'Preferences');
  % make sure we do not got into infinite loop/crash
  fig = mifit_fig;
  set(fig, 'CloseRequestFcn','closereq');  % sets to default
  
  % save the window size and location
  set(fig,'Units','pixels');
  config.Position = get(fig,'Position');
  mifit_Preferences_Save(config);
  
  if isfield(config, 'Save_Data_On_Exit') && strcmp(config.Save_Data_On_Exit, 'yes')
    mifit_File_Save;
  else
    file = fullfile(prefdir, [ mfilename '.mat' ]);
    if ~isempty(dir(file)), delete(file); end
  end
  mifit_disp([ '[Exit] Exiting miFit. Bye bye.' ]);
  if ishandle(mifit_fig('mifit_View_Parameters'))
    delete(mifit_fig('mifit_View_Parameters'))
  end
  % clearing the Models appdata may lead to a crash
  rmappdata(fig, 'Models');
  rmappdata(fig, 'DnDControl');
  close(fig);
  
