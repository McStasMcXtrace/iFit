function mifit_File_Saveas(varargin)
% File/saveas: save as a MAT with user filename
% save the application configuration into specified file
% Data sets and Model parameters into a mifit.mat file
  fig        = mifit_fig;
  Data       = getappdata(fig, 'Data');
  Models     = getappdata(fig, 'Models');
  Optimizers = getappdata(fig, 'Optimizers');
  CurrentOptimizer         = getappdata(fig, 'CurrentOptimizer');
  CurrentOptimizerConfig   = getappdata(fig, 'CurrentOptimizerConfig');
  CurrentOptimizerCriteria = getappdata(fig, 'CurrentOptimizerCriteria');
  CurrentModel             = getappdata(fig, 'CurrentModel');
 
  if nargin == 1 && ischar(varargin{1})
    file = varargin{1};
  else
    filterspec = { '*.mat','MAT-files (*.mat)'};
    [filename, pathname] = uiputfile(filterspec, 'Save All miFit Data sets and Models as', [ 'mifit.mat' ]);
    if isequal(filename,0) || isequal(pathname,0)
      return
    end
    file = fullfile(pathname, filename);
  end
  builtin('save', file, 'Data','Models','Optimizers','CurrentOptimizer', ...
    'CurrentOptimizerConfig', 'CurrentOptimizerCriteria');
  d = dir(file);
  mifit_disp([ '[File_Saveas] Saved Data sets/Models into ' file ' [size: ' num2str(round(d.bytes/1024)) ' kb]' ]);
