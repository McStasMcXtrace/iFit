function [pars,criteria,message,output] = fits(a, model, pars, options, constraints,ub)
% [pars,criteria,message,output] = fits(a, model, pars, options, constraints) : fit data set on a model
%
%   @iData/fits find best parameters estimates in order to minimize the 
%     fitting criteria using function 'fun' as model, by mean of an optimization method
%     described with the 'options' structure argument.
%     Additional constraints may be set by fxing some parameters, or define
%     more advanced constraints (min, max, steps). The last arguments controls the fitting
%     options with the optimset mechanism, and the constraints to apply during optimization.
%  [pars,...] = fits(a, model, pars, options, lb, ub)
%     uses lower and upper bounds as parameter constraints (double arrays)
%  [pars,...] = fits(a, model, pars, options, fixed)
%     indicates which parameters are fixed (non zero elements of array).
%  [pars,...] = fits(a, model, pars, 'optimizer', ...)
%     uses a specific optimizer and its default options.
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
% options.optimizer
%   Optimization method. Default is 'fminsearch' (char/function handle)
%   the syntax for calling the optimizer is e.g. optimizer(criteria,pars,options,constraints)
% options.criteria
%   Minimization criteria. Default is 'least_square' (char/function handle)
%   the syntax for evaluating the criteria is criteria(Signal, Error, Model)
% options.OutputFcn
%   Function called at each iteration as outfun(pars, optimValues, state)
%   The 'fminplot' function may be used.
% options.Display
%   Display additional information during fit: 'iter','off','final'. Default is 'iter'.
%
% input:  a: object or array (iData)
%         model: model function (char/cellstr)
%         pars: initial model parameters (double array)
%         options: structure as defined by optimset/optimget (char/struct)
%           if given as a char, it defines the algorithm to use and its default options
%         constraints: fixed parameter array. Use 1 for fixed parameters, 0 otherwise (double array)
%           OR use a structure with some of the following fields:
%           constraints.min:   minimum parameter values (double array)
%           constraints.max:   maximum parameter values (double array)
%           constraints.step:  maximum parameter step/change allowed.
%           constraints.fixed: fixed parameter flag. Use 1 for fixed parameters, 0 otherwise (double array)
%
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
%         o=optimset('fminsearch'); o.OutputFcn='fminplot'; 
%         [p,c,m,o]=fits(a,'gauss',[1 2 3 4],o);
%
% Version: $Revision: 1.16 $
% See also iData, fminsearch, optimset, optimget

% nested  functions: outfun_wrapper, eval_criteria
% private functions: least_square, fits_constraints
 
% handle default parameters, if missing
if nargin < 3
  pars = [];
end
if nargin < 5
  constraints = [];
end
if nargin < 4, options=[]; end
if isempty(options)
  options = 'fminsearch';
end
if ischar(options) | isa(options, 'function_handle')
  algo = options;
  options           = feval(algo,'defaults');
  options.optimizer = algo;
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
if ~isfield(options, 'optimizer')
  options.optimizer = 'fminsearch';
end
if ~isfield(options, 'criteria')
  options.criteria  = @least_square;
end

[dummy, info] = ieval(a, model,'identify'); % model info 
if isempty(pars)
  pars=info.Guess;               % get default starting parameters
end
if isnumeric(constraints) | islogical(constraints)
  if nargin<6
    fixed            = constraints;
    constraints      =[];
    constraints.fixed=fixed;
  else
    lb = constraints;
    constraints.min = lb;
    constraints.max = ub;
  end
end
if ~isstruct(constraints)
  iData_private_error(mfilename,[ 'The constraints argument is of class ' class(constraints) '. Should be a single array or a struct' ]);
end

% removes warnings
try
  warn.set = warning('off','iData:setaxis');
  warn.get = warning('off','iData:getaxis');
catch
  warn = warning('off');
end
% remove output function calls so that we use our own, and overload with the user's choice
if isfield(options, 'OutputFcn'); options.UserOutputFcn = options.OutputFcn; 
else options.UserOutputFcn=''; end
options.OutputFcn     = @outfun_wrapper;

if ~isfield(options,'Display') options.Display=''; end
if ~isfield(options,'algorithm') options.algorithm=options.optimizer; end
if isempty(options.Display)    options.Display='notify'; end

pars = reshape(pars, [ 1 numel(pars)]); % a single row
constraints.parsStart      = pars;
constraints.parsHistory    = [];
constraints.parsNames      = info.Parameters;
constraints.criteriaHistory= [];
constraints.modelName      = info.Name;
constraints.algorithm      = options.algorithm;
constraints.optimizer      = options.optimizer;
constraints.funcCounts     = 0;

if strcmp(options.Display, 'iter') | strcmp(options.Display, 'final') | strcmp(options.Display, 'notify')
  disp([ '** Starting fit of ' a.Tag ' using model ' info.Name ' with optimizer ' options.algorithm ]);
  disp(char(a))
  disp(  '** Minimization performed on parameters:');
  disp(info.Parameters(:)');
end

t0 = clock;

% call minimizer
try
[pars_out,criteria,message,output] = feval(options.optimizer, ...
    @(pars) eval_criteria(pars, model, options.criteria, a), pars, options, constraints);
catch
[pars_out,criteria,message,output] = feval(options.optimizer, ...
    @(pars) eval_criteria(pars, model, options.criteria, a), pars, options);
end

% apply constraints on final results, if any
pars = fits_constraints(pars_out, constraints);

% set output/results
if ischar(message) | ~isfield(output, 'message')
  output.message = message;
else
  output.message = [ '(' num2str(message) ') ' output.message ];
end
output.optimizer  = options.optimizer;
output.funcCounts = constraints.funcCounts;
output.modelName  = constraints.modelName;
output.modelInfo  = info;
output.modelValue = ieval(a, model, pars);
output.pars       = pars;
output.parsNames  = constraints.parsNames;
output.parsHistory= constraints.parsHistory;
output.criteriaHistory=constraints.criteriaHistory;
output.criteria   = criteria;
output.duration   = etime(clock, t0);
output.options    = options;
output.constraints= constraints;


% reset warnings
try
  warning(warn.set);
  warning(warn.get);
catch
  warning(warn);
end

% Use a nested function as the OutputFcn wrapper, to access 'constraints', 'options' and 'output'
  function stop = outfun_wrapper(pars, optimValues, state);
    % we need to transform pars first
    pars = fits_constraints(pars, constraints);
    stop = false;
    
    % then call the user supplied OutputFcn
    if ~isempty(options.UserOutputFcn)
      stop1 = feval(options.UserOutputFcn, pars, optimValues, state);
      stop = stop | stop1;
    end
        
    if isfield(options, 'MaxFunEvals')
    if options.MaxFunEvals > 0 & constraints.funcCounts >= 1.2*options.MaxFunEvals
      stop=true;
    end
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
    m      = get(a,'Monitor'); m=real(m);
    if not(all(m == 1) | all(m == 0)),
      Signal = Signal./m; Error=Error./m; % per monitor
    end
    
    % compute criteria
    c = feval(criteria, Signal(:), Error(:), Model(:));
    
    % save current optimization state
    constraints.criteriaHistory = [ constraints.criteriaHistory ; c ];
    constraints.parsPrevious    = pars;
    constraints.parsHistory     = [ constraints.parsHistory ; pars ]; 
    constraints.funcCounts      = constraints.funcCounts+1; 
  end
  
end % fits end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c=least_square(Signal, Error, Model)
% weighted least square criteria, which is also the Chi square
  if all(Error == 0)
    c = sum(abs(Signal-Model).^2); % raw least square
  else
    index = find(Error & ~isnan(Error) & ~isinf(Error));
    c=(Signal(index)-Model(index))./Error(index);
    c=abs(c);
    c=sum(c.*c);                % Chi square
  end
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
pars=pars(:)';
end
