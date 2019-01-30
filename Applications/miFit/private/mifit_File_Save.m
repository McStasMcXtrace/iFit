function mifit_File_Save(varargin)
% File/Save: save Data sets and Model parameters into a mifit.mat file
  
  % save Data and Models
  file = fullfile(prefdir, [ 'mifit.mat' ]);
  mifit_File_Saveas(file);
  
  % save Preferences with window size and location
  config = getappdata(mifit_fig, 'Preferences');
  set(fig,'Units','pixels');
  config.Position = get(fig,'Position');
  mifit_Preferences_Save(config);
