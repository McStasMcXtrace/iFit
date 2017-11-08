function mifit_Data_Fit(varargin)
% Data/Fit: fit selected data sets with their attached Model
  [d, index_selected] = mifit_List_Data_pull;
  if isempty(index_selected) || numel(d) == 0, return; end
  id = [ d{:} ]; % a single iData array
  if all(isempty(id)), return; end
  
  % get the Optimizer configuration
  CurrentOptimizer = getappdata(mifit_fig,'CurrentOptimizer');
  if isappdata(mifit_fig,'CurrentOptimizerConfig')
    options = getappdata(mifit_fig,'CurrentOptimizerConfig');
  else
    options = [];
  end
  % overload Preferences choices: OutputFcn, Display, Criteria.
  config = getappdata(mifit_fig, 'Preferences');
  options.Display='iter';
  if isfield(config,'Fit_Verbose') && ...
          (strcmp(config.Fit_Verbose,'yes') || strcmp(config.Fit_Verbose,'on') ...
          || strcmp(config.Fit_Verbose,'1') || ...
            (isscalar(config.Fit_Verbose) && config.Fit_Verbose == 1))
    if ~isfield(options, 'PlotFcns') || isempty(options.PlotFcns), options.PlotFcns = {}; end
    options.PlotFcns{end+1} = 'fminplot';
    options.PlotFcns{end+1} = @(x,optimValues,state)mifit_Models_View_Parameters(x);
  else
    options.OutputFcn = 'fminstop';
  end
  if isappdata(mifit_fig,'CurrentOptimizerCriteria')
    options.criteria=getappdata(mifit_fig,'CurrentOptimizerCriteria');
  end
  if isfield(options, 'criteria') && isempty(options.criteria)
    options = rmfield(options, 'criteria');
  end
  
  set(mifit_fig,'Pointer','watch');
  
  % [pars,criteria,message,output] = fits(a, model, pars, options, constraints, ...)
  mifit_disp([ 'Starting fit of data set(s) with "' CurrentOptimizer '"' ]);
  mifit_disp(char(id))
  
  % ********* THE FIT *********
  % the initial Dataset array 'd' is updated after the fit.
  try
    [p,c,m,o]=fits(id, '', 'current', options);  % with assigned models or gaussians
  catch ME
    disp(getReport(ME))
    set(mifit_fig,'Pointer','arrow');
    return
  end

  % update Data list with fit results (and History)
  if numel(id) == 1
    d = { id };
  else
    for index=1:numel(id)
      d{index} = id(index);
    end
  end
  
  D = getappdata(mifit_fig, 'Data');
  if numel(D) == 1
    D = d;
  else
    D(index_selected) = d;
  end
  setappdata(mifit_fig, 'Data',D);
  mifit_History_push;
  
  % log fit results to the miFit log file
  for index=1:numel(d)
    
    this_d = d{index}; 
    if numel(d) == 1
      output = o; this_p     = p;
    else 
      output = o{index}; this_p     = p{index};
    end
    mifit_disp([ '** Final fit results for ' char(this_d) ], true);
    sigma = output.parsHistoryUncertainty;
    if ~isempty(output.parsHessianUncertainty)
      jj=output.constraints.index_variable;
      sigma(jj) = max(sigma(jj), output.parsHessianUncertainty);
    end
    if isfield(output, 'constraints')
      constraints= output.constraints;
    else constraints = []; end

    for indexp=1:numel(p)
      t0=sprintf('  p(%3d):%20s=%g +/- %g', indexp,strtok(output.parsNames{indexp}), this_p(indexp), sigma(indexp)); 
      if isfield(constraints, 'fixed') && length(constraints.fixed) >= indexp && constraints.fixed(indexp)
        t1=' (fixed)'; else t1=''; end
      mifit_disp([ t0 t1 ], true);  % only in the Log file.
    end
    if isfield(output,'corrcoef')
      mifit_disp(sprintf(' Correlation coefficient=%g (closer to 1 is better)',  output.corrcoef), true);
      mifit_disp(sprintf(' Weighted     R-factor  =%g (Rwp, smaller that 0.2 is better)', output.Rfactor), true);
    end
    
    if isfield(output, 'parsHessianCorrelation') && ~isempty(output.parsHessianCorrelation)
      corr = output.parsHessianCorrelation;
      nb_true_independent_parameters = sum(1./sum(corr.^2));
      mifit_disp([ ' Estimated number of independent parameters: ' num2str(nb_true_independent_parameters) ], true)
      mifit_disp(' Correlation matrix (non diagonal terms indicate non-independent parameters):', true)
      mifit_disp(corr, true);
    end
    
    % update plot
    % evaluate model with its parameters (Edit) and Data set axes
    model = this_d.Model;
    if model.Duration < 0.5
      % update plot, if found
      h = mifit_fig([ 'plot_' this_d.Tag ]);
      if ~isempty(h)
        plot(this_d,'light transparent grid tight replace');
      end
    end
  end
  
  % show results for the last fit/dataset
  index_selected = index_selected(end); 
  
  setappdata(mifit_fig, 'CurrentDataSetIndex', index_selected);
  setappdata(mifit_fig, 'CurrentDataSet', this_d);
  
  % update Parameter Window content
  handle = mifit_Models_View_Parameters('update');
  if ~isempty(handle) && ishandle(handle)
    setappdata(handle, 'LastFitOutput', o);
  end
  
  % display the Parameter distribution histograms
  if isfield(config,'Fit_Verbose') && strcmp(config.Fit_Verbose,'yes')
    mifit_Models_View_Parameters('histograms');
  end
  set(mifit_fig,'Pointer','arrow');
  
