function [pars,criteria,message,output] = fits(a, model, pars, constraints, options)
% [pars,criteria,message,output] = fits(a, model, pars, constraints, options) : fit data set on a model
%
%   @iData/fits find best parameters estimates in order to minimize the 
%     fitting criteria using function 'fun' as model, by mean of an optimization method.
%     Additional constraints may be set by fxing some parameters, or define
%     more advanced constraints (min, max, steps). The last argument controls the fitting
%     options with the optimset mechanism.
%
% options.TolX
%   The termination tolerance for x. Its default value is 1.e-4.
% options.TolFun
%   The termination tolerance for the function value. The default value is 1.e-4. 
%   This parameter is used by fminsearch, but not fminbnd.
% options.MaxIter
%   Maximum number of iterations allowed.
% options.MaxFunEvals
%   The maximum number of function evaluations allowed. 
%   The default value is 500 and 200*length(x0) for fminsearch
% options.algorithm
%   Optimization method. Default is 'fminsearch' (char/function handle)
%   the syntax for calling the optimizer is e.g. algorithm(criteria,pars,options)
% options.criteria
%   Minimization criteria. Default is 'least_square' (char/function handle)
%   the syntax fo evaluating the criteria is criteria(Signal, Error, Model)
% options.OutputFcn
%   Function called at each iteration as outfun(pars, optimValues, state)
% options.PlotFcns
%   Similar to OutputFcn, but dedicated to plotting. 
%   The 'fits_plot' function may be used.
% options.Display
%   Display additional information during fit 'iter','off','final'
%
% input:  a: object or array (iData)
%         model: model function (char/cellstr)
%         pars: initial model parameters (double array)
%         constraints: fixed parameter array. Use 1 for fixed parameters, 0 otherwise (double array)
%           OR use a structure with some of the following fields:
%           constraints.min:   minimum parameter values (double array)
%           constraints.max:   maximum parameter values (double array)
%           constraints.step:  maximum parameter step/change allowed.
%           constraints.fixed: fixed parameter flag. Use 1 for fixed parameters, 0 otherwise (double array)
%         options: structure as defined by optimset/optimget (char/struct)
% output: 
%         pars:              best parameter estimates (double array)
%         criteria:          minimal criteria value achieved (double)
%         message:           return message/exitcode from the optimizer (char)
%         output:            additional information about the optimization (structure)
%           algorithm:         Algorithm used (char)
%           funcCount:         Number of function evaluations (double)
%           iterations:        Number of iterations (double)
%           parsHistory:       Parameter set history during optimization (double array)
%           criteriaHistory:   Criteria history during optimization (double array)
%           modelValue:        Last model evaluation (iData)
%
% ex:     p=fits(a,'gauss',[1 2 3 4]);
%         o=optimset('fminsearch'); o.PlotFcns='fits_plot'; 
%         [p,c,m,o]=fits(a,'gauss',[1 2 3 4],[],o);
%
% See also iData, fminsearch, optimset, optimget
 
% handle default parameters, if missing
if nargin < 3
  pars = [];
end
if nargin < 4
  constraints = [];
end
if nargin < 5
  options = optimset('fzero');
end

% handle input iData arrays
if length(a) > 1
 pars_out={}; criteria=[]; message={}; output={};
  for index=1:length(a(:))
    [pars_out{index}, criteria(index), message{index}, output{index}] = ...
      fits(a(index), model, pars, constraints, options);
  end
  pars = pars_out;
  return
end

% handle options
if ~isfield(options, 'algorithm')
  options.algorithm = @fminsearch;
end
if ~isfield(options, 'criteria')
  options.criteria  = @least_square;
end

[dummy, info] = ieval(iData, model); % model info 
if isempty(pars)
  pars=info.Guess;               % get default starting parameters
end
if isnumeric(constraints) | islogical(constraints)
  tmp              = double(constraints);
  constraints      =[];
  constraints.fixed=tmp;
end
if ~isstruct(constraints)
  iData_private_error(mfilename,[ 'The constraints argument is of class ' class(constraints) '. Should be a single array or a struct' ]);
end
 
constraints.parsStart=pars;

% removes warnings
try
  warn.set = warning('off','iData:setaxis');
  warn.get = warning('off','iData:getaxis');
catch
  warn = warning('off');
end
% remove output function calls so that we use own own, and overload with the user's choice
options.UserOutputFcn = options.OutputFcn;
options.OutputFcn     = @outfun_wrapper;
if isfield(options, 'PlotFcns')
  options.UserPlotFcns= options.PlotFcns;
  options.PlotFcns    = [];
end
pars = reshape(pars, [ 1 numel(pars)]); % a single row
constraints.parsHistory    = [];
constraints.criteriaHistory= [];
constraints.parsNames      = info.Parameters;
constraints.modelName      = info.Name;
constraints.algorithm      = options.algorithm;

% call minimizer
[pars_out,criteria,message,output] = feval(options.algorithm, ...
    @(pars) eval_criteria(pars, model, options.criteria, a), pars, options);

% apply constraints, if any
pars = fits_constraints(pars_out, constraints);

% set output/results
output.criteriaHistory=constraints.criteriaHistory;
output.parsHistory= constraints.parsHistory;
output.modelValue = ieval(a, model, pars);
output.criteria   = criteria;
outout.pars       = pars;
output.message    = [ '(' num2str(message) ') ' output.message ];
output.algorithm  = options.algorithm;

% reset warnings
try
  warning(warn.set);
  warning(warn.get);
catch
  warning(warn);
end

% Use a nested function as the OutputFcn wrapper, to access 'constraints'
  function stop = outfun_wrapper(pars, optimValues, state);
    % we need to transform pars first
    pars = fits_constraints(pars, constraints);
    stop = false;
    
    % then call the user supplied OutputFcn
    if ~isempty(options.UserOutputFcn)
      stop1 = feval(options.UserOutputFcn, pars, optimValues, state);
      stop = stop | stop1;
    end
    
    % call the optional Plotting function
    if ~isempty(options.UserPlotFcns)
      optimValues.parsHistory    = constraints.parsHistory;
      optimValues.criteriaHistory= constraints.criteriaHistory;
      optimValues.parsNames      = constraints.parsNames;
      optimValues.algorithm      = constraints.algorithm;
      optimValues.modelName      = constraints.modelName;
      try
        stop2 = feval(options.UserPlotFcns, pars, optimValues, state);
      catch
        lasterr
        stop2=true;
      end
      stop = stop | stop2;
    end
    
  end

% Use a nested function as the criteria wrapper, to access 'constraints'
  function c = eval_criteria(pars, model, criteria, a)
  % criteria to minimize
  
    % apply constraints on pars first
    pars = fits_constraints(pars,constraints);
    
    % then get data
    Signal = get(a,'Signal');
    Error  = get(a,'Error');
    Model  = ieval(a, model, pars);
    Model  = get(Model, 'Signal');
    
    % compute criteria
    c = feval(criteria, Signal(:), Error(:), Model(:));
    
    % save current optimization state
    constraints.criteriaHistory = [ constraints.criteriaHistory ; c ];
    constraints.parsPrevious    = pars;
    constraints.parsHistory     = [ constraints.parsHistory ; pars ]; 
    
  end
  
end % fits end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c=least_square(Signal, Error, Model)
% weighted least square criteria, which is also the Chi square
if all(Error == 0)
  c = sum((Signal-Model).^2); % raw least square
else
  index = find(Error);
  c=(Signal(index)-Model(index))./Error(index);
  c=abs(c);
  c=sum(c.*c);                % Chi square
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stop = fits_plot(pars, optimValues, state)
% default plotting function showing the criteria evolution 
% as well as main parameters and status

% fields in optimValues: funcount, fval, iteration, procedure, parsHistory, criteriaHistory
% state values: 'init','interrupt','iter','done'
  stop = false;
  h = findall(0, 'Tag', 'iData_fits_plot');
  if length(h) > 1, delete(h(2:end)); h=h(1); end
  if isempty(h) & optimValues.iteration <=1
    h = figure('Tag','iData_fits_plot', 'Unit','pixels');
    tmp = get(h, 'Position'); tmp(3:4) = [500 400];
    set(h, 'Position', tmp);
  elseif isempty(h)
    stop = true;  % figure was closed: abort optimization by user
  end
  figure(h);
  set(h, 'Name', [ 'iData/fits: ' func2str(optimValues.algorithm) ' iteration ' num2str(optimValues.iteration) ]);
  subplot(1,2,1); % this subplot shows the criteria
  plot(optimValues.criteriaHistory);
  xlabel('iteration'); ylabel('criteria'); axis tight
  if strcmp(state, 'done'), title('Done'); 
  elseif strcmp(state, 'init'), title('Init'); 
  else title('Close figure to abort');  end
  
  subplot(1,2,2); % this subplot shows some parameters
  pars_hist= optimValues.parsHistory;
  pars_std = std(pars_hist)./mean(pars_hist);
  [dummy, sort_std] = sort(pars_std);
  pars_hist = pars_hist(:, sort_std);
  switch length(pars_std)
  case 1
    plot(pars_hist, optimValues.criteriaHistory);
    xlabel([ optimValues.parsNames{1} ]);
    ylabel('criteria');
  case 2
    fscatter3(pars_hist(:,1), pars_hist(:,2), optimValues.criteriaHistory);
    xlabel([ num2str(sort_std(1)) ': ' optimValues.parsNames{sort_std(1)} ]);
    ylabel([ num2str(sort_std(2)) ': ' optimValues.parsNames{sort_std(2)} ]);
    zlabel('criteria'); view(3);
  otherwise
    fscatter3(pars_hist(:,1), pars_hist(:,2), pars_hist(:,3), optimValues.criteriaHistory);
    xlabel([ num2str(sort_std(1)) ': ' optimValues.parsNames{sort_std(1)} ]);
    ylabel([ num2str(sort_std(2)) ': ' optimValues.parsNames{sort_std(2)} ]);
    zlabel([ num2str(sort_std(3)) ': ' optimValues.parsNames{sort_std(3)} ]);
    view(3); c = colorbar;
  end
  title({[ func2str(optimValues.algorithm) ' iteration ' num2str(optimValues.iteration) ],optimValues.modelName});
  axis tight
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pars = fits_constraints(pars, constraints)
% take into account constraints on fit parameters

if isfield(constraints, 'fixed')  % fix some parameters
  index = find(constraints.fixed);
  pars(index) = constraints.parsStart(index);
else
  if isfield(constraints, 'min')  % lower bound for parameters
    index = find(pars < constraints.min);
    pars(index) = constraints.min(index);
  end
  if isfield(constraints, 'max')  % upper bound for parameters
    index = find(pars > constraints.max);
    pars(index) = constraints.max(index);
  end
  if isfield(constraints, 'step') % restrict parameter change
    parsStep = pars - constraints.parsPrevious;
    index = find( abs(parsStep) > abs(constraints.steps) );
    parsStep(index) = sign(parsStep).*abs(constraints.steps);
    pars(index) = constraints.parsPrevious(index) + parsStep(index);
  end
end
end
