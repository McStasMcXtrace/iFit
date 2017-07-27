function [pars,fval,exitflag,output] = fmin_private_wrapper(optimizer, fun, pars, options, constraints, varargin)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fmin_private_wrapper(OPTIMIZER, FUN,PARS,[OPTIONS],[constraints], ...) wrapper to optimizers
%
%  Checks for input arguments and options. Then calls the optimizer with a wrapped 
%  inline_objective function, which applies constraints and makes stop condition checks.
%  the main optimizer call is within a try/catch block which exits when an early 
%  stop is met.
% 
% Calling:
%   fmin_private_wrapper(optimizer, fun, pars) asks to minimize the 'fun' inline_objective function with starting
%     parameters 'pars' (vector)
%   fmin_private_wrapper(optimizer, fun, pars, options) same as above, with customized options (optimset)
%   fmin_private_wrapper(optimizer, fun, pars, options, fixed) 
%     is used to fix some of the parameters. The 'fixed' vector is then 0 for
%     free parameters, and 1 otherwise.
%   fmin_private_wrapper(optimizer, fun, pars, options, lb, ub) 
%     is used to set the minimal and maximal parameter bounds, as vectors.
%   fmin_private_wrapper(optimizer, fun, pars, options, constraints) 
%     where constraints is a structure (see below).
%   fmin_private_wrapper(optimizer, problem) where problem is a structure with fields
%     problem.inline_objective:   function to minimize
%     problem.x0:          starting parameter values
%     problem.options:     optimizer options (see below)
%     problem.constraints: optimization constraints
%   fmin_private_wrapper(optimizer, ..., args, ...)
%     sends additional arguments to the objective function
%       criteria = FUN(pars, args, ...)
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fmin_private_wrapper('fminimfil',banana,[-1.2, 1])
%
% Input:
%  OPTIMIZER is the name/handle to an optimizer function, or '' for default
%
%  FUN is the function to minimize (handle or string): criteria = FUN(PARS)
%  It needs to return a single value or vector.
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess. PARS can also be given as a single-level structure.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%     o=fmin_private_wrapper(optimizer,'defaults')
%  An empty OPTIONS sets the default configuration.
%
%  constraints may be specified as a structure
%     constraints.min=   vector of minimal values for parameters
%     constraints.max=   vector of maximal values for parameters
%     constraints.fixed= vector having 0 where parameters are free, 1 otherwise
%     constraints.step=  vector of maximal parameter changes per iteration
%  An empty CONSTRAINTS sets no constraints.
%
%  Additional arguments are sent to the objective function.
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% Version: $Date$
% See also: fminsearch, optimset

% NOTE: all optimizers have been gathered here so that maintenance is minimized
% each user call function only defines the options... The optimizer by itself is
% in the 'private'.
%
% private: 'inline_objective', 
%          'inline_apply_constraints', 
%          'inline_constraints_minmax', 
%          'inline_private_check'
%          'inline_localChar', 
%          'inline_disp'
%          'inline_estimate_uncertainty'
%          'inline_auto_optimizer'

% return code     message
%  0                Algorithm terminated normally
% -1                Termination function tolerance criteria reached
% -2                Maximum number of iterations reached
% -3                Maximum number of function evaluations reached
% -4                Function value is Inf or Nan
% -5                Termination parameter tolerance criteria reached
% -6                Algorithm was terminated by the output function
% -7                Maximum consecutive rejections exceeded (anneal)
% -8                Minimum temperature reached (anneal)
% -9                Global Simplex convergence reached (simplex)
% -10               Optimization terminated: Stall Flights Limit reached (swarm)
% -11               Other termination status (cmaes/Ralg)
% -12               Termination function change tolerance criteria reached

% input parameter handling =====================================================

% nargin stuff (number of parameters)
% default options for optimset

if nargin < 1,         optimizer =''; end
if nargin < 2,         fun       = ''; end
if nargin < 3,         pars      = []; end
if isempty(optimizer), optimizer = 'fmin'; end
if isempty(fun),       fun       = 'defaults'; end

if nargin < 4
  options=[];
end
if nargin < 5
  constraints = [];
end
if nargin < 6
  varargin = {};
end
objective = fun;
if strcmp(fun,'defaults') || strcmp(fun,'identify')
  pars = feval(optimizer, 'defaults');
  return
elseif nargin < 2
  error([ 'syntax is: ' optimizer '(optimizer, objective, parameters, ...)' ] );
elseif nargin >= 2 && isstruct(fun)
  if     isfield(fun, 'x0'),          pars=fun.x0;
  elseif isfield(fun, 'guess'),       pars=fun.guess;
  elseif isfield(fun, 'Guess'),       pars=fun.Guess; end
  if     isfield(fun, 'options'),     options=fun.options; end
  if     isfield(fun, 'constraints'), constraints=fun.constraints; end 
  if     isfield(fun, 'objective'),   tmp=fun.objective; fun=[]; fun=tmp; 
  elseif isfield(fun, 'model'),       tmp=fun.model; fun=[]; fun=tmp;
  elseif isfield(fun, 'f'),           tmp=fun.f; fun=[]; fun=tmp;
  elseif isfield(fun, 'function'),    tmp=fun.function; fun=[]; fun=tmp; end
elseif nargin < 3
  error([ 'syntax is: ' inline_localChar(optimizer) '(objective, parameters, ...)' ] );
end
if isempty(pars)
  error([ inline_localChar(optimizer) ': starting parameters (3rd argument) must not be empty.' ] );
end

if ~ischar(fun) && ~isa(fun, 'function_handle')
  error([ inline_localChar(optimizer) ': objective function (2nd argument) must be a char or function_handle, but is a ' class(fun) '.' ] );
end


% default arguments when missing ===============================================
if isempty(options)
  options=optimizer;
end

if (ischar(options) && exist(options) == 2) | isa(options, 'function_handle')
  optimizer = options;
  options=feval(options, 'defaults');
  options.optimizer = optimizer;
elseif ischar(options), options=str2struct(options); end
if ischar(pars),
  pars   =str2struct(pars); 
end

% handle case when parameters are given as structures
% constraints =============================================================

if length(constraints)==length(pars) & (isnumeric(constraints) | islogical(constraints))
  if nargin < 6,               % given as fixed index vector
    fixed             = constraints; 
    constraints       = [];
    constraints.fixed = fixed;  % avoid warning for variable redefinition.
  elseif isnumeric(varargin{1}) && ~isempty(varargin{1}) ...
      && length(constraints) == length(varargin{1})
    % given as lb,ub parameters (nargin==6)
    lb = constraints; 
    ub = varargin{1};
    varargin(1) = []; % remove the 'ub' from the additional arguments list
    constraints     = [];
    constraints.min = lb;
    constraints.max = ub;
  end
end

if isstruct(pars)
  pars_name=fieldnames(pars);
  pars=cell2mat(struct2cell(pars));
  % check if constraints are also structures, but with same fields
  check = {'min','max','fixed','step'};
  for index=1:length(check)
    if isfield(constraints,check{index}) && isstruct(constraints.(check{index}))
      new = [];
      for f=1:length(pars_name)
        if isfield(constraints.(check{index}), pars_name{f})
          new = [ new constraints.(check{index}).(pars_name{f}) ];
        end
      end
      if length(new) == length(pars)
        constraints.(check{index}) = new;
      else
        error([ inline_localChar(optimizer) ': parameters and constraints %s are given as structures, but not with same length/fields (%i and %i resp.).' ],  ...
          check{index}, length(pars_name), length(fieldnames(constraints.(check{index}))) );
      end
    end
  end
else
  pars_name = {};
end

if length(constraints)==length(pars) & (isnumeric(constraints) | islogical(constraints))
  if nargin < 6,               % given as fixed index vector
    fixed             = constraints; 
    constraints       = [];
    constraints.fixed = fixed;  % avoid warning for variable redefinition.
  elseif isnumeric(varargin{1}) && ~isempty(varargin{1}) ...
      && length(constraints) == length(varargin{1})
    % given as lb,ub parameters (nargin==6)
    lb = constraints; 
    ub = varargin{1};
    varargin(1) = []; % remove the 'ub' from the additional arguments list
    constraints     = [];
    constraints.min = lb;
    constraints.max = ub;
  end
end

if ~isempty(constraints) && ischar(constraints)
  constraints = str2struct(constraints);
end
if ~isempty(constraints) && ~isstruct(constraints)
  error([ inline_localChar(optimizer) ': The constraints argument is of class ' class(constraints) '. Should be a vector or a struct' ]);
end
if ~isstruct(options)
  error([ inline_localChar(optimizer) ': The options argument is of class ' class(options) '. Should be a string or a struct' ]);
end

pars = pars(:)';
constraints.parsStart       = pars;  % used when applying constraints
if ~isfield(constraints,'fixed')
  constraints.fixed=zeros(size(pars));
end
if ~isfield(constraints,'pars_name')
  constraints.pars_name=pars_name;
end

% reduce the number of free parameters to only those which are not fixed
pars_all = pars;
constraints.index_variable = find(~constraints.fixed | isnan(constraints.fixed));
pars = pars(constraints.index_variable);
% NOTE: from here 'pars' only contains the variable parameters

constraints.parsPrevious    = pars;
constraints.parsBest        = pars;
constraints.parsHistory     = [];
constraints.criteriaHistory = [];
constraints.criteriaStart   = [];
constraints.criteriaPrevious= Inf;
constraints.criteriaBest    = Inf;
constraints.funcCount       = 0;
constraints.message         = '';
constraints.fevalDuration   = 0;
options.pars_name           = constraints.pars_name;  % avail for inline_disp

if isfield(constraints, 'Expression') && ~isfield(constraints, 'eval')
  constraints.eval = constraints.Expression;
end

% options ======================================================================
if ~isfield(options,'optimizer') || isempty(options.optimizer)
  options.optimizer = optimizer;
end
if ~isfield(options,'Display') options.Display=''; end

options=inline_private_check(options, feval(options.optimizer,'defaults'));
t0=clock;

n = prod(size(pars)); N=n; % these are for the following 'eval's
numberOfVariables = n;
numberofvariables = n;
if ischar(options.MaxFunEvals), 
  options.MaxFunEvals = eval(options.MaxFunEvals); 
end
if ischar(options.MinFunEvals), 
  options.MinFunEvals = eval(options.MinFunEvals); 
end
if ischar(options.MaxIter), 
  options.MaxIter = eval(options.MaxIter); 
end

if ischar(options.TolFun)
  options.TolFunChar = options.TolFun;
  if options.TolFun(end)=='%'
    options.TolFun(end)='';
    fval = inline_objective_all(fun, pars_all, options, constraints, varargin{:});
    options.TolFun = abs(str2num(options.TolFun)*fval/100);
  else
    options.TolFun = str2num(options.TolFun);
  end
else
  fval = NaN;
end

if ischar(options.TolX)
  options.TolXChar=options.TolX;
  if options.TolX(end)=='%'
    options.TolX(end)='';
    options.TolX = abs(str2num(options.TolX)*pars(:)/100);
  else
    options.TolX = str2num(options.TolX);
  end
end
if ~isfield(options, 'algorithm'), options.algorithm = options.optimizer; end

% Optimizer call ===============================================================
if strncmp(options.Display,'iter',4)
  disp([ '** Starting minimization of ' inline_localChar(fun) ' using algorithm ' inline_localChar(options.algorithm) ]);
  if ~isempty(constraints.pars_name)
    spars = sprintf(' %s', constraints.pars_name{constraints.index_variable});
  else spars = ''; end
  disp([ 'Func_count  min[f(x)]    Parameters' spars ]);
  inline_disp(options, constraints.funcCount , pars, fval)
end

message    = constraints.message;
exitflag   = 0;
iterations = 0;
fval       = Inf;       % in case this is a vector, it should be a row
pars       = pars(:);   % should be a column for most optimizers
output     = [];

try

% calls the optimizer with a wrapped 'inline_objective' function
%    which applies constraints and makes stop condition checks.
% the main optimizer call is within a try/catch block which exits when an early 
%  stop is met. See private 'inline_objective' and 'inline_apply_constraints'
%
  inline_call_optimizer(fun, pars_all, pars, options, constraints, varargin{:});
  
catch ME
  % we may get here when the inline_objective issue an error after a stop condition.
  output.lasterror = ME;
  message = output.lasterror.message;
end % try

if isstruct(output) && isfield(output,'lasterror') && isempty(strfind(output.lasterror.message, 'stop condition:'))
  disp('Code error when launching the optimizer. Please fix it...')
  disp(output.lasterror.message);
  for index=1:length(output.lasterror.stack)
    disp(output.lasterror.stack(index))
  end
  % this is a real error (not from a stop condition)
  rethrow(output.lasterror);
end

% post optimization checks =====================================================

fval = constraints.criteriaBest; % set in inline_objective
fval = sum(fval(:));
pars = constraints.parsBest;
pars = pars(:)';  % returned as a row

if iterations, 
  output.iterations    = iterations;
elseif ~isfield(output,'iterations')
  output.iterations    = constraints.funcCount ;
end

% determine message (return status of optimiser)
if isempty(message) && isfield(constraints,'message')
  message = constraints.message;
end

if isempty(message)
  if exitflag==0
    message='Algorithm terminated';
  end
end
if ~isfield(output,'message')
  output.message         = message;
end

output.funcCount       = constraints.funcCount ;
output.algorithm       = options.algorithm;
output.parsHistory     = constraints.parsHistory;
output.criteriaHistory = constraints.criteriaHistory;
output.parsBest        = constraints.parsBest;
output.criteriaBest    = constraints.criteriaBest;
output.options         = options; 
output.constraints     = constraints;
output.optimizer       = options.optimizer;
output.duration        = etime(clock, t0);
output.fevalDuration   = constraints.fevalDuration;
output.parsStart       = constraints.parsStart;

% estimate parameter uncertainty from the search trajectory ====================

index      = find(output.criteriaHistory < min(output.criteriaHistory)*4);   % identify tolerance region around optimum 
if length(index) < 3 % retain 1/4 lower criteria part
  delta_criteria = output.criteriaHistory - min(output.criteriaHistory);
  index      = find(abs(delta_criteria/min(output.criteriaHistory)) < 0.25);
end
if length(index) < 3
  index = 1:length(output.criteriaHistory);
end
try
  delta_pars = (output.parsHistory(index,:)-repmat(output.parsBest,[length(index) 1])); % get the corresponding parameter set
  weight_pars= exp(-((output.criteriaHistory(index)-min(output.criteriaHistory))).^2 / 8); % Gaussian weighting for the parameter set
  weight_pars= repmat(weight_pars,[1 length(output.parsBest)]);
  output.parsHistoryUncertainty = sqrt(sum(delta_pars.*delta_pars.*weight_pars)./sum(weight_pars));
end

covp = [];
if ((strcmp(options.Display,'final') || strcmp(options.Display,'iter') ...
  || (strcmp(options.Display,'notify') && isempty(strfind(message, 'Converged')))) || nargout == 4) ...
  && ((isfield(options,'Diagnostics') && strcmp(options.Diagnostics,'on')) ...
 || (length(pars)^2*output.fevalDuration < 60 ... % should spend less than a minute to compute the Hessian
  && (~isfield(options,'Diagnostics') || ~strcmp(options.Diagnostics,'off')) ...
  && exitflag ~= -6)) % not when user explicitely requested premature end (Abort)
  if length(pars)^2*output.fevalDuration > 5
    disp([ '  Estimating Hessian matrix... (' num2str(length(pars)^2*output.fevalDuration) ' [s] remaining, please wait)' ]);
  end
  try
    [dp, covp, corp,jac,hessian]  = ...
      inline_estimate_uncertainty(fun, pars, options, varargin{:});
    if ~isempty(covp)
      output.parsHessianUncertainty = reshape(abs(dp), size(pars));
      output.parsHessianCovariance  = covp;
      output.parsHessianCorrelation = corp;
      output.parsHessian            = hessian;
      output.parsJacobian           = jac;
    end
  end
end
if isempty(covp)
  output.parsHessianUncertainty = [];
  output.parsHessianCovariance  = [];
  output.parsHessianCorrelation = [];
  output.parsHessian            = [];
  output.parsJacobian           = [];
end

if strcmp(options.Display,'final') || strcmp(options.Display,'iter') ...
  || (strcmp(options.Display,'notify') && isempty(strfind(message, 'Converged'))) ...
  || (isfield(options,'Diagnostics') && strcmp(options.Diagnostics,'on'))
  disp([ sprintf('\n') '** Finishing minimization of ' inline_localChar(fun) ' using algorithm ' inline_localChar(options.algorithm) ]);
  disp( [ ' Status: ' output.message ]);
  if ~isempty(constraints.pars_name)
    spars = sprintf(' %s', constraints.pars_name{constraints.index_variable});
  else spars = ''; end
  disp([ ' Func_count     min[f(x)]        Parameters' spars ]);
  inline_disp(struct('Display','iter'), -constraints.funcCount , pars, mean(fval));
  
  % test length of tolerance region
  if length(index) > 10 && isfield(output,'parsHistoryUncertainty') 
    disp(' Gaussian uncertainty on parameters (half width, from the optimization history)')
    inline_disp(struct('Display','iter'), -1, output.parsHistoryUncertainty, NaN);
  end
  if isfield(output,'parsHessianUncertainty') && ~isempty(output.parsHessianUncertainty)
    disp(' Gaussian uncertainty on parameters (half width, from the Hessian matrix)')
    inline_disp(struct('Display','iter'), -1, output.parsHessianUncertainty, NaN);
  end
  if isfield(output, 'parsHessianCorrelation') && ~isempty(output.parsHessianCorrelation)
    % the trace of the correlation matrix is M. The non diagonal terms indicate 
    % correlated parameters
    corr = output.parsHessianCorrelation;
    n_free = numel(pars);
    if isfield(constraints, 'fixed')  % fix some parameters
      n_free = n_free - numel(find(constraints.fixed & ~isnan(constraints.fixed)));
      corr=corr(find(~constraints.fixed),find(~constraints.fixed));
    end
    nb_true_independent_parameters = sum(1./sum(corr.^2));
    
    if n_free > ceil(nb_true_independent_parameters*1.2)+1
      warn = [ '. WARNING: too many parameters (' num2str(numel(n_free)) ') in problem' ]; 
    else warn = ''; end
    disp([ ' Estimated number of independent parameters: ' num2str(nb_true_independent_parameters) warn ])
    disp(' Correlation matrix (non diagonal terms indicate non-independent parameters):')
    disp(corr)
  end
end

% restore initial parameters as a structure (when given as such)
if ~isempty(constraints.pars_name)
  pars_name = constraints.pars_name(:);
  pars = cell2struct(num2cell(pars(:)), ...
    strtok(pars_name(constraints.index_variable)), 1);
else % provide full parameter vector
  pall = constraints.parsStart;
  
  index_variable = find(~constraints.fixed | isnan(constraints.fixed));
  pall(index_variable) = pars;
  pars = pall;
end

% ==============================================================================
%     |_   _| \ | | |    |_   _| \ | |  ____| 
%       | | |  \| | |      | | |  \| | |__    
%       | | | . ` | |      | | | . ` |  __|   
%      _| |_| |\  | |____ _| |_| |\  | |____  
%     |_____|_| \_|______|_____|_| \_|______| 
% ==============================================================================
% inline function which carry the constraints, options, pars, fval, etc
% even when an error is triggered, e.g. when a stop condition is raised.

  function inline_call_optimizer(fun, pars_all, pars, options, constraints, varargin)
% This is where we call the actual optimizer. 
% This is controlled from fmin_private_wrapper

  if strcmp(options.optimizer, 'fmin')
    [optimizer, algorithm] = inline_auto_optimizer(fun, pars_all, varargin{:});
    % update options
    options2 =feval(optimizer, 'defaults');
    % merge the structures/cells handling duplicated fields
    options = mergestruct(options2, options);
    options.optimizer = optimizer;
    options.algorithm = algorithm;
  end
  
  switch options.optimizer
  case {'cmaes','fmincmaes'}    
  % Evolution Strategy with Covariance Matrix Adaption ---------------------------
    hoptions.MaxIter    = options.MaxIter;
    hoptions.TolFun     = options.TolFun;
    hoptions.MaxFunEvals= options.MaxFunEvals;
    hoptions.PopSize    = options.PopulationSize;
    hoptions.DispFinal  = 'off';
    hoptions.DispModulo = 0;
    hoptions.SaveVariables  = 'off';
    hoptions.LogModulo      = 0;
    hoptions.LogPlot        = 'off';
    if isfield(constraints,'step'), hoptions.DiffMaxChange = constraints.step(:); end
    if isfield(constraints,'min'),  hoptions.LBounds=constraints.min(:); end
    if isfield(constraints,'max'),  hoptions.UBounds=constraints.max(:); end
    if isfield(constraints,'min') & isfield(constraints,'max')
      sigma = abs(constraints.max(:) - constraints.min(:))/4;
    elseif isfield(constraints,'step')
      sigma = constraints.step*10;
    else
      sigma = 0.3;
    end

    [pars, fval, iterations, exitflag, output] = cmaes(@(pars) inline_objective(fun, pars, varargin{:}), pars, ...
      sigma, hoptions);

    if     strmatch(exitflag, 'tolx')
      exitflag=-5;
      message = [ 'Converged: Termination parameter tolerance criteria reached (options.TolX=' ...
                num2str(options.TolX) ')' ];
    elseif strmatch(exitflag, 'tolfun')
      exitflag=-1;
      message = [ 'Converged: Termination function tolerance criteria reached (options.TolFun=' ...
                num2str(options.TolFun) ')' ];
    elseif strmatch(exitflag, 'maxiter')
      exitflag=-2;
      message = [ 'Maximum number of iterations reached (options.MaxIter=' ...
                num2str(options.MaxIter) ')' ];
    elseif strmatch(exitflag, 'maxfunevals')
      exitflag=-3;
      message = [ 'Maximum number of function evaluations reached (options.MaxFunEvals=' ...
                num2str(options.MaxFunEvals) ')' ];
    elseif strmatch(exitflag, 'outputfcn')
      exitflag=-6;
      message = 'Algorithm was terminated by the output function (options.OutputFcn)';
    elseif strmatch(exitflag, 'funvalcheck')
      exitflag=-4;
      message = 'Function value is Inf or Nan (options.FunValCheck)';  
    else   
      message=char(exitflag);
      exitflag=-11;
    end
  case {'ga','fminga','GA'}          
  % genetic algorithm ------------------------------------------------------------
    [constraints, constraints_var] = inline_constraints_minmax(pars_all, constraints);
    [pars,fval,iterations,output] = GA(@(pars) inline_objective(fun, pars, varargin{:}), pars, options,constraints_var);
  case {'gradrand','ossrs','fmingradrand'}
  % random gradient --------------------------------------------------------------
    [pars,fval,iterations] = ossrs(pars, @(pars) inline_objective(fun, pars, varargin{:}), options);
  case {'hooke','fminhooke'}
  % Hooke-Jeeves direct search ---------------------------------------------------
    [pars,histout] = hooke(pars, @(pars) inline_objective(fun, pars, varargin{:}), ...
                         options.MaxFunEvals, 2.^(-(0:options.MaxIter)), options.TolFun);
    iterations      = size(histout,1);
  case {'imfil','fminimfil'}
  % Unconstrained Implicit filtering (version 1998) ------------------------------
    [pars,fval,iterations,output] = imfil(pars, @(pars) inline_objective(fun, pars, varargin{:}), options); 
  case {'fminlm','LMFsolve'}
  % Levenberg-Maquardt steepest descent ------------------------------------------
    % LMFsolve minimizes the sum of the squares of the inline_objective: sum(inline_objective.^2)
    [pars, fval, iterations, exitflag] = LMFsolve(@(pars) inline_objective(fun, pars, varargin{:}), pars, ...
             'Display',0, 'FunTol', options.TolFun, 'XTol', options.TolX, ...
             'MaxIter', options.MaxIter, 'Evals',options.MaxFunEvals);
    switch exitflag
    case -5, message='Converged: Termination function tolerance criteria reached';
    case -2, message='Maximum number of iterations reached';
    case -3, message='Maximum number of function evaluations reached';
    case -1, message='Converged: Termination parameter tolerance criteria reached';
    end
  case {'powell','fminpowell'}
  % Powell minimization ----------------------------------------------------------
    if isempty(options.Hybrid), options.Hybrid='Coggins'; end
    if strcmp(lower(options.Hybrid), 'coggins') 
      t = 'Coggins';
    else 
      t = 'Golden rule';
    end
    options.algorithm  = [ 'Powell Search (by Secchi) [' options.optimizer '/' t ']' ];
    [constraints, constraints_var] = inline_constraints_minmax(pars_all, constraints);
    [pars,fval,exitflag,output] = powell(@(pars) inline_objective(fun, pars, varargin{:}), pars, options);
  case {'pso','fminpso'}
  % particle swarm ---------------------------------------------------------------
    [constraints, constraints_var] = inline_constraints_minmax(pars_all, constraints);
    [pars,fval,exitflag,output] = PSO(@(pars) inline_objective(fun, pars, varargin{:}),pars, ...
       constraints_var.min(:),constraints_var.max(:),options);
    message = output.message;
  case {'ralg','fminralg','solvopt'}
  % Shor's r-algorithm -----------------------------------------------------------
    opt(1) = -1;
    opt(2) = options.TolX;
    opt(3) = options.TolFun;
    opt(4) = options.MaxIter;
    opt(5) = -1;
    opt(6) = 1e-8; 
    opt(7) = 2.5; 
    opt(8) = 1e-11;

    % call the optimizer
    [pars,fval,out,iterations, message] = ralg(pars, @(pars) inline_objective(fun, pars, varargin{:}), ...
      [], opt,[],[], [], options.MaxFunEvals, options.FunValCheck);
    if out(9) < 0, exitflag = out(9); 
    else exitflag=0; end
  case {'fminsearch','fminsearchbnd'}
  % Nelder-Mead simplex, with constraints ----------------------------------------
    [pars,fval,exitflag,output] = fminsearch(@(pars) inline_objective(fun, pars, varargin{:}), pars, options);
  %     1  Maximum coordinate difference between current best point and other
  %        points in simplex is less than or equal to TolX, and corresponding 
  %        difference in function values is less than or equal to TolFun.
  %     0  Maximum number of function evaluations or iterations reached.
  %    -1  Algorithm terminated by the output function.
    if     exitflag == 1, exitflag=-5;
    elseif exitflag == 0, exitflag=-3;
    elseif exitflag ==-1, exitflag=-6;
    end
  case {'simpsa','fminsimpsa','SIMPSA'}
  % simplex/simulated annealing --------------------------------------------------
    [constraints, constraints_var] = inline_constraints_minmax(pars_all, constraints);
    [pars,fval,exitflag,output] = SIMPSA(@(pars) inline_objective(fun, pars, varargin{:}), pars, ...
      constraints_var.min(:),constraints_var.max(:),options);
  case {'SCE','fminsce'}
  % shuffled complex evolution ---------------------------------------------------
    [constraints, constraints_var] = inline_constraints_minmax(pars_all, constraints);
    [pars,fval,exitflag,output] = SCE(@(pars) inline_objective(fun, pars, varargin{:}), pars, ...
      constraints_var.min(:),constraints_var.max(:),options);
  case {'hPSO','fminswarmhybrid','fminswarm'}
    [constraints, constraints_var] = inline_constraints_minmax(pars_all, constraints);
    if isa(options.Hybrid, 'function_handle') | exist(options.Hybrid) == 2
      hoptions.algorithm = [ 'hybrid Particle Swarm Optimizer (by Leontitsis) [' options.optimizer '/' inline_localChar(options.Hybrid) ']' ];
    else
      hoptions.algorithm = [ 'Particle Swarm Optimizer (by Leontitsis) [fminswarm]' ];
    end

    % transfer optimset options and constraints
    hoptions.space     = [ constraints_var.min(:) constraints_var.max(:) ];
    hoptions.MaxIter   = options.MaxIter;
    hoptions.TolFun    = options.TolFun;
    hoptions.TolX      = options.TolX;
    hoptions.Display   = options.Display;
    hoptions.MaxFunEvals=options.MaxFunEvals;
    hoptions.FunValCheck=options.FunValCheck;
    hoptions.OutputFcn  =options.OutputFcn;

    hoptions.Hybrid     =options.Hybrid;
    hoptions.c1         =options.SwarmC1;
    hoptions.c2         =options.SwarmC2;
    hoptions.w          =options.SwarmW;
    hoptions.bees       =options.PopulationSize;
    if isa(options.Hybrid, 'function_handle') | exist(options.Hybrid) == 2
      hoptions.StallFliLimit = 50;
    else
      hoptions.StallFliLimit = Inf;
    end
    if isfield(constraints,'step')
      hoptions.maxv = constraints.step;
    else
      hoptions.maxv = abs(constraints_var.max(:)-constraints_var.min(:))/2;
    end
    [pars,fval,iterations,output] = hPSO(@(pars) inline_objective(fun, pars, varargin{:}), pars, hoptions);
  case {'Simplex','fminsimplex'}
  % Nelder-Mead simplex state machine --------------------------------------------
    [constraints, constraints_var] = inline_constraints_minmax(pars_all, constraints);
    [pars, out]=Simplex('init', pars, abs(constraints_var.max(:)-constraints_var.min(:))/10);  % Initialization
    for iterations=1:options.MaxIter
      fval = feval(@(pars) inline_objective(fun, pars, varargin{:}), pars);

      [pars,out]=Simplex( fval );
      if isfield(options,'TolFunChar')
        options.TolFun = options.TolFunChar;
        if options.TolFun(end)=='%'
          options.TolFun(end)='';
          options.TolFun = abs(str2num(options.TolFun)*fval/100);
        else
          options.TolFun = str2num(options.TolFun);
        end
      end
      if Simplex('converged', options.TolFun)             % Test for convergence
        exitflag=-1;
        message= [ 'Converged: Termination function tolerance criteria reached (options.TolFun=' ...
                  num2str(options.TolFun) ')' ];
        break
      end
      if iterations >= options.MaxIter
        exitflag=-2;
        message = [ 'Maximum number of iterations reached (options.MaxIter=' ...
                  num2str(options.MaxIter) ')' ];;
        break
      end
    end
    pars=Simplex('centroid'); % obtain the final value.
  case {'cgtrust','fmincgtrust'}
  % Steihaug Newton-CG-Trust region algorithm ------------------------------------
    [pars,histout] = cgtrust(pars, @(pars) inline_objective(fun, pars, varargin{:}), ...
      [ options.TolFun .1 options.MaxIter options.MaxIter], options.TolX*options.TolX);
  % [pars,histout] = levmar(pars, @(pars) inline_objective(fun, pars, varargin{:}), options.TolFun, options.MaxIter);
    iterations      = size(histout,1);
  % not so efficient optimizers ==================================================
  case {'fminanneal','anneal'}  
  % simulated annealing ----------------------------------------------------------
    options.MaxTries   = options.MaxIter/10;
    options.StopVal    = options.TolFun;
    options.Verbosity=0;
    [pars,fval,iterations,exitflag] = anneal(@(pars) inline_objective(fun, pars, varargin{:}), pars(:)', options);
    if exitflag==-7,message='Maximum consecutive rejections exceeded (anneal)'; end
  case {'fminbfgs','bfgs'}      
  % Broyden-Fletcher-Goldfarb-Shanno ---------------------------------------------
    [pars, histout, costdata,iterations] = bfgswopt(pars(:), @(pars) inline_objective(fun, pars, varargin{:}), options.TolFun, options.MaxIter);
    iterations = size(histout,1);
  case {'fminkalman','kalmann','ukfopt'}
  % unscented Kalman filter ------------------------------------------------------
    [pars,iterations] = ukfopt(@(pars) inline_objective(fun, pars(:)), pars(:), ...
                options.TolFun, norm(pars)*eye(length(pars)), 1e-6*eye(length(pars)), 1e-6);
  case {'ntrust','fminnewton'}
  % Dogleg trust region, Newton model --------------------------------------------
    [pars,histout,costdata] = ntrust(pars(:),@(pars) inline_objective(fun, pars, varargin{:}), ...
         options.TolFun,options.MaxIter);
    iterations      = size(histout,1);
  case {'buscarnd','fminrand'}
  % adaptive random search -------------------------------------------------------
    [pars,fval]=buscarnd(@(pars) inline_objective(fun, pars, varargin{:}), pars, options);
  case {'markov','mcmc','gwmcmc','fminmarkov'}
    % make sure parameters are bound
    [constraints, constraints_var] = inline_constraints_minmax(pars_all, constraints);
    % create initial distribution
    M=numel(pars); %number of model parameters
    Nwalkers=max(2*M, 40); %number of walkers/chains.
    minit=randn(M,Nwalkers);
    LB=constraints_var.min;
    UB=constraints_var.max;
    for index=1:M
      minit(index,:) = LB(index) + rand(1,Nwalkers).*(UB(index)-LB(index));
    end
    % set nb of iterations
    mccount = max(options.MaxFunEvals, options.MaxIter);
    % clean options
    if isfield(options,'StepSize')
      op.StepSize = options.StepSize; else op.StepSize = 2; end
    if isfield(options,'ThinChain')
      op.ThinChain = options.ThinChain; else op.ThinChain = 10; end
    if isfield(options,'ProgressBar')
      op.ProgressBar = options.ProgressBar; else op.ProgressBar = true; end
    if isfield(options,'BurnIn')
      op.BurnIn = options.BurnIn; else op.BurnIn = 0; end
    % launch the MCMC. The criteria is assumed as a log, negative and to be maximized
    [pars,fval]=gwmcmc(minit, @(pars) -inline_objective(fun, pars, varargin{:}), ...
      mccount, op);
    pars = pars(:,:);
    fval = fval(:,:);
    
  otherwise
    % unknown optimizer. 
    error([ inline_localChar(optimizer) ': Unknown optimizer.' ]);
    return
  end % switch
  
  % we may never get there if a stop condition triggers and error
  % this error is then recovered by fmin_private_wrapper in the try catch around
  % the optimizer call.
  
  end % inline_call_optimizer




  % ============================================================================
  % inline objective function so that we can use options and constraints
  function fval = inline_objective(fun, p, varargin)
    % prepares an objective function which only has variable parameters
  
    % new objective fun(p, varargin)
    pall = constraints.parsStart;
    p    = p(:)';
    
    index_variable = find(~constraints.fixed | isnan(constraints.fixed));
    pall(index_variable) = p;
    
    t = clock;
    % call the objective with all parameters
    [fval, pall] = inline_objective_all(fun, pall, options, constraints, varargin{:});
    p = pall(index_variable);
    
    % save current optimization state
    if sum(fval) < sum(constraints.criteriaBest(:)), 
      constraints.criteriaBest=fval;
      constraints.parsBest    =p;
    end
    
    constraints.parsPrevious    = pall;
    constraints.parsHistory     = [ constraints.parsHistory ; p(:)' ]; 
    
    % check for usual stop conditions MaxFunEvals, TolX, TolFun ..., and call OutputFcn
    [exitflag, message] = inline_private_check(pall, sum(fval), ...
       constraints.funcCount , options, constraints);
    constraints.message = message;
    
    if isempty(constraints.criteriaStart)
      constraints.criteriaStart = fval;
    end
    constraints.fevalDuration   = etime(clock, t); % time required to estimate the criteria
    constraints.criteriaPrevious= fval;
    constraints.criteriaHistory = [ constraints.criteriaHistory ; sum(constraints.criteriaPrevious)
     ];
    constraints.funcCount       = constraints.funcCount +1; 
    
    if exitflag
      error([ 'stop condition: ' message ]); % will end optimization in try/catch
    end
    
    inline_disp(options,  constraints.funcCount, p, fval);
  end % inline_objective

end % fmin_private_wrapper
