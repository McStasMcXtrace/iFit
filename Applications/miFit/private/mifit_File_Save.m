function mifit_File_Save(varargin)
% File/Save: save Data sets and Model parameters into a mifit.mat file
  
  file = fullfile(prefdir, [ 'mifit.mat' ]);
  mifit_File_Saveas(file);
