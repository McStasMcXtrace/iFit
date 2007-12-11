function [pars,criteria,message,output] = fits(a, model, pars, constraints, options)
% [pars,criteria,message,output] = fits(a, model, pars, constraints, options) : fit data set on a model
%
%   @iData/fits find best parameters estimates in order to minimize the 
%     fitting criteria using function 'fun' as model, by mean of an optimization method.
%
% input:  a: object or array (iData)
%         model: model function (char/cellstr)
%         pars: initial model parameters (double array)
%         constraints: fixed parameter array. Use 1 for fixed parameters, 0 otherwise (double array or structure)
%           constraints.min: minimum parameter values (double array)
%           constraints.max: maximum parameter values (double array)
%           constraints.step:maximum parameter step allowed. Use 0 for fixed parameters (double array)
%           constraints.fixed:fixed parameter flag. Use 1 for fixed parameters, 0 otherwise (double array)
%         options: structure as defined by optimset/optimget (char/struct)
%           options.algorithm: optimization method. Default is 'fminsearch' (char/function handle)
%             the syntax for calling the optimizer is e.g. fminsearch(criteria,pars,options)
%           options.criteria: minimization criteria. Default is 'least_square' (char/function handle)
%             the syntax fo evaluating the criteria is least_square(Signal, Error, Model)
% output: information about the optimization (structure)
%         pars:              best parameter estimates (double array)
%         criteria:          minimal criteria value achieved (double)
%         message:           return message/exitcode from the optimizer (char)
%         algorithm:         Algorithm used (char)
%         funcCount:         Number of function evaluations (double)
%         iterations:        Number of iterations (double)
%         parameterHistory:  Parameter set history during optimization (double array)
%         criteriaHistory:   Criteria history during optimization (double array)
%         modelValue:        Last model evaluation (iData)
% ex:     p=fits(a,'gauss',[1 2 3 4]);
%
% See also iData, fminsearch, optimset, optimget

% exitflag returns a string that describes the exit condition of optimizer
%     OK: Optimizer converged to a solution.
%     OK: Algorithm was terminated by the output function.
%     OK: Change in x was smaller than the specified tolerance.
%     OK: Change in the objective function value was less than the specified tolerance.
%     OK: Magnitude of gradient smaller than the specified tolerance.
%     WARNING: Line search cannot find an acceptable point along the current search direction.
%     WARNING: Maximum number of function evaluations or iterations was reached.
%     WARNING: Number of iterations exceeded options.MaxIter
%     WARNING: Number of function evaluations exceeded options.FunEvals.
%     ERROR: Bounds are inconsistent (x1 > x2).
%
%options.Display
%  A flag that determines if intermediate steps in the minimization appear on the screen. 
%  If set to 'iter', intermediate steps are displayed; if set to 'off', no intermediate 
%  solutions are displayed, if set to final, displays just the final output.
%options.TolX
%  The termination tolerance for x. Its default value is 1.e-4.
%options.TolFun
%  The termination tolerance for the function value. The default value is 1.e-4. 
%  This parameter is used by fminsearch, but not fminbnd.
%options.MaxIter
%  Maximum number of iterations allowed.
%options.MaxFunEvals
%  The maximum number of function evaluations allowed. The default value is 500 for 
%  fminbnd and 200*length(x0) for fminsearch
%options.OutputFcn
%  Specify a user-defined function that the optimization function calls at each iteration.
%options.FunValCheck
%  Check whether objective function values are valid. 'on' displays a warning when 
%  the objective function returns a value that is complex, Inf or NaN. 'off' (the default) displays no warning.
%options.PlotFcns         
%  Plots various measures of progress while the algorithm
%  executes, select from predefined plots or write your
%  own. Specifying @optimplotx plots the current point;
%  @optimplotfunccount plots the function count;
%  @optimplotfval plots the function value;
%  @optimplotstepsize plots the step size;
%  @optimplotfirstorderopt plots the first-order of
%  optimality.
%options.Diagnostics
%  Display diagnostic information about the function to be minimized.
%options.DiffMaxChange
%  Maximum change in variables for finite differencing.
%options.DiffMinChange
%  Minimum change in variables for finite differencing.
 
% handle default parameters, if missing
% handle input options: depends on chosen method
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

if ~isfield(options, 'algorithm')
  options.algorithm = @fminsearch;
end
if ~isfield(options, 'criteria')
  options.criteria = @least_square;
end
if isempty(pars)
  [dummy, pars] = ieval(a, model);
end

% removes warnings
try
  warn.set = warning('off','iData:setaxis');
  warn.get = warning('off','iData:getaxis');
catch
  warn = warning('off');
end

% call minimizer
[pars_out,criteria,message,output] = feval(options.algorithm, ...
    @(pars) eval_criteria(pars, model, options.criteria, a), pars, options);

% set output/results
pars = pars_out;
output.modelValue = ieval(a, model, pars);
output.criteria   = criteria;
outout.pars       = pars;
output.message    = message;
output.algorithm  = options.algorithm;

% reset warnings
try
  warning(warn.set);
  warning(warn.get);
catch
  warning(warn);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = eval_criteria(pars, model, criteria, a)
% criteria to minimize

  Signal = get(a,'Signal');
  Error  = get(a,'Error');
  Model  = ieval(a, model, pars);
  Model  = get(Model, 'Signal');
  
  % return criteria
  c = feval(criteria, Signal(:), Error(:), Model(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c=least_square(Signal, Error, Model)
% weighted least square criteria, which is also the Chi square
if all(Error == 0)
  c = sum((Signal-Model).^2);
else
  index = find(Error);
  c=(Signal(index)-Model(index))./Error(index);
  c=abs(c);
  c=sum(c.*c);
end

