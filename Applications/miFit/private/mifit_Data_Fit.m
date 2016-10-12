function mifit_Data_Fit(varargin)
% Data/Fit: fit selected data sets with their attached Model
  [d, index_selected] = mifit_List_Data_pull;
  if isempty(index_selected), return; end
  
  % get the Optimizer configuration
  CurrentOptimizer = getappdata(mifit_fig,'CurrentOptimizer');
  if isappdata(mifit_fig,'CurrentOptimizerConfig')
    options = getappdata(mifit_fig,'CurrentOptimizerConfig');
  else
    options = [];
  end
  % overload Preferences choices: OutputFcn, Display.
  config = getappdata(mifit_fig, 'Preferences');
  options.Display='iter';
  if isfield(config,'Fit_Verbose') && strcmp(config.Fit_Verbose,'yes')
    if ~isfield(options, 'PlotFcns') || isempty(options.PlotFcns), options.PlotFcns = {}; end
    options.PlotFcns{end+1} = 'fminplot';
    options.PlotFcns{end+1} = @(x,optimValues,state)mifit_Models_View_Parameters(x);
  end
  
  % [pars,criteria,message,output] = fits(a, model, pars, options, constraints, ...)
  mifit_disp([ 'Starting fit of data set(s) with "' CurrentOptimizer '"' ]);
  mifit_disp(char(d))
  
  % ********* THE FIT *********
  % the initial Dataset array 'd' is updated after the fit.
  [p,c,m,o]=fits(d, '', 'current', options);  % with assigned models or gaussians

  % update Data list with fit results (and History)
  D = getappdata(mifit_fig, 'Data');
  if numel(D) == 1
    D = d;
  else
    D(index_selected) = d;
  end
  setappdata(mifit_fig, 'Data',D);
  mifit_History_push;
  
  % show results for the 1st fit/dataset
  index_selected = index_selected(1); 
  d=d(1);
  
  setappdata(mifit_fig, 'CurrentDataSetIndex', index_selected);
  setappdata(mifit_fig, 'CurrentDataSet', d);
  
  % update Parameter Window content
  handle = mifit_Models_View_Parameters('update');
  if ~isempty(handle) && ishandle(handle)
    setappdata(handle, 'LastFitOutput', o);
  end
  
  % display the Parameter distribution histograms
  if isfield(config,'Fit_Verbose') && strcmp(config.Fit_Verbose,'yes')
    mifit_Models_View_Parameters('histograms');
  end
