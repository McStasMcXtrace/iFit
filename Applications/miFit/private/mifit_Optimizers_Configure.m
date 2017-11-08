function mifit_Optimizers_Configure(varargin)
% Optimizers/Configure: open dialogue to change current optimizer parameters
  % change optimizer configuration parameters
  CurrentOptimizer = getappdata(mifit_fig,'CurrentOptimizer');
  mifit_disp([ '[Optimizers_Configure] Configure "' CurrentOptimizer '"' ]);
  options    = getappdata(mifit_fig, 'CurrentOptimizerConfig');
  config     = getappdata(mifit_fig, 'Preferences');
  o.FontSize = config.FontSize;
  o.Name     = [ mfilename ': Optimizer "' CurrentOptimizer '" configuration' ];
  NL = sprintf('\n');
  o.TooltipString = [ 'Most relevant "' CurrentOptimizer '" configuration items:' NL ...
    '* MaxFunEvals - Maximum number of function evaluations allowed [ positive integer ]' NL ...
    '* MaxIter - Maximum number of iterations allowed [ positive scalar ]' NL ...
    '* TolFun - Termination tolerance on the function value [ positive scalar ]' NL ...
    '* TolX - Termination tolerance on X [ positive scalar ]' ];
  options1 = uitable(options, o); % structure GUI as a Table (spreadsheet)
  if isempty(options1), return; end
  % look for changes in new options...
  fields = fieldnames(options);
  for index=1:numel(fields)
    if ~isequal(options1.(fields{index}), options.(fields{index}))
      t1 = class2str(fields{index},options1.(fields{index})); t1 = strrep(t1, sprintf('\n'), '');
      t0 = class2str(fields{index},options.(fields{index}));  t0 = strrep(t0, sprintf('\n'), '');
      mifit_disp([ '[Optimizers_Configure] Assigned ' CurrentOptimizer ': ' ...
        t1 ' [was ' t0 ']' ]);
    end
  end
  options = options1;
  setappdata(mifit_fig,'CurrentOptimizerConfig', options);
