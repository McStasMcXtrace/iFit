function [pars,fval,exitflag,output] = fmin_private_wrapper(optimizer, fun, pars, options, constraints, ub)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fmin_private_wrapper(OPTIMIZER, FUN,PARS,[OPTIONS],[CONSTRAINTS]) wrapper to optimizers
%
%  Checks for input arguments and options. Then calls the optimizer with a wrapped 
%  objective function, which applies constraints and makes stop condition checks.
%  the main optimizer call is within a try/catch block which exits when an early 
%  stop is met.
% 
% Calling:
%   fmin_private_wrapper(optimizer, fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fmin_private_wrapper(optimizer, fun, pars, options) same as above, with customized options (optimset)
%   fmin_private_wrapper(optimizer, fun, pars, options, fixed) 
%     is used to fix some of the parameters. The 'fixed' vector is then 0 for
%     free parameters, and 1 otherwise.
%   fmin_private_wrapper(optimizer, fun, pars, options, lb, ub) 
%     is used to set the minimal and maximal parameter bounds, as vectors.
%   fmin_private_wrapper(optimizer, fun, pars, options, constraints) 
%     where constraints is a structure (see below).
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fmin_private_wrapper('fminimfil',banana,[-1.2, 1])
%
% Input:
%  OPTIMIZER is the name/handle to an optimizer function, or '' for default
%  FUN is a function handle (anonymous function or inline) with a loss
%  function, which may be of any type, and needn't be continuous. It does,
%  however, need to return a single value.
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the simulated annealing, 
%  compliant with optimset. Default options may be obtained with
%     o=fminanneal('defaults')
%  options.TEMP_START sets the starting temperature
%  options.TEMP_END   sets the end temperature
%
%  CONSTRAINTS may be specified as a structure
%   constraints.min=   vector of minimal values for parameters
%   constraints.max=   vector of maximal values for parameters
%   constraints.fixed= vector having 0 where parameters are free, 1 otherwise
%   constraints.step=  vector of maximal parameter changes per iteration
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% Version: $Revision: 1.1 $
% See also: fminsearch, optimset

% NOTE: all optimizers have been gathered here so that maintenance is minimized
% each user call function only defined the options...
%
% private: 'objective', 'apply_constraints', 'constraints_minmax'

% parameter handling ===========================================================

% nargin stuff (number of parameters)
% default options for optimset

if nargin < 1, optimizer=''; end
if nargin < 2, fun = ''; end
if isempty(optimizer), optimizer = 'fminimfil'; end
if isempty(fun),       fun = 'defaults'; end

if strcmp(fun,'defaults')
  pars = feval(optimizer, 'defaults');
  return
elseif nargin < 2
  error([ 'syntax is: ' mfilename '(optimizer, objective, parameters, ...)' ] );
elseif nargin < 3
  error([ 'syntax is: ' localChar(optimizer) '(objective, parameters, ...)' ] );
end
if nargin < 4
  options=[];
end
if nargin < 5
  constraints = [];
end
if nargin < 6
  ub = [];
end

% default arguments when missing
if isempty(options)
  options=feval(optimizer, 'defaults');
end
if length(constraints) & isnumeric(constraints) 
  if nargin == 4,               % given as fixed index vector
    fixed = constraints; constraints=[];
    constraints.fixed = fixed;  % avoid warning for variable redefinition.
  else                          % given as lb,ub parameters (nargin==5)
    lb = constraints; clear constraints;
    constraints.min = lb;
    constraints.max = ub;
  end
end

if isfield(constraints, 'min')  % test if min values are valid
  index=find(isnan(constraints.min));
  constraints.min(index) = -2*abs(pars(index));
end
if isfield(constraints, 'max')  % test if max values are valid
  index=find(isnan(constraints.max));
  constraints.max(index) = 2*abs(pars(index));
end

constraints.parsStart       = pars;  % used when applying constraints
constraints.parsPrevious    = pars;
constraints.parsBest        = pars;
constraints.parsHistory     = [];
constraints.criteriaHistory = [];
constraints.criteriaStart   = [];
constraints.criteriaPrevious= Inf;
constraints.criteriaBest    = Inf;
constraints.funcCounts      = 0;

options=fmin_private_std_check(options, feval(options.optimizer,'defaults'));
t0=clock;

n = prod(size(pars));
numberOfVariables = n;
numberofvariables = n;
if ischar(options.MaxFunEvals), 
  options.MaxFunEvals = eval(options.MaxFunEvals); 
end
if ischar(options.MaxIter), 
  options.MaxIter = eval(options.MaxIter); 
end

if strcmp(options.Display,'iter')
  fmin_private_disp_start(mfilename, fun, pars);
end

message    = 'Algorithm terminated normally';
exitflag   = 0;
iterations = 0;
fval=Inf;
output=[];

% Optimizer call ===============================================================

try
% calls the optimizer with a wrapped 'objective' function
%    which applies constraints and makes stop condition checks.
% the main optimizer call is within a try/catch block which exits when an early 
%  stop is met. See private 'objective' and 'apply_constraints' below.

switch options.optimizer
case {'fminanneal','anneal'}  
% simulated annealing ----------------------------------------------------------
  options.MaxTries   = options.MaxIter;
  options.StopVal    = options.TolFun;
  switch options.Display
  case 'iter',  options.Verbosity=2;
  case 'final', options.Verbosity=1;
  otherwise,    options.Verbosity=0;
  end
  [pars,fval,iterations] = anneal(@(pars) objective(fun, pars), pars, options);
case {'fminbfgs','bfgs'}      
% Broyden-Fletcher-Goldfarb-Shanno ---------------------------------------------
  [pars, histout, costdata,iterations] = bfgswopt(pars(:), @(pars) objective(fun, pars), options.TolFun, options.MaxIter);
  fval = constraints.criteriaBest;
  iterations = size(histout,1);
case {'cmaes','fmincmaes'}    
% Evolution Strategy with Covariance Matrix Adaption ---------------------------
  hoptions.MaxIter    = options.MaxIter;
  hoptions.TolFun     = options.TolFun;
  hoptions.MaxFunEvals= options.MaxFunEvals;
  %hoptions.FunValCheck= options.FunValCheck;
  %hoptions.OutputFcn  = options.OutputFcn;
  hoptions.PopSize    = options.PopulationSize;
  %hoptions.algorithm  = options.algorithm;
  if strcmp(options.Display,'final'), hoptions.DispFinal = 'on';
  else                                hoptions.DispFinal = 'off'; end
  if strcmp(options.Display,'iter'),  hoptions.DispModulo=1; 
  else hoptions.DispModulo=0; end
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
  hoptions.SaveVariables  = 'off';
  hoptions.LogModulo      = 0;
  hoptions.LogPlot        = 'off';
  
  [pars, fval, iterations, exitflag, output] = cmaes(@(pars) objective(fun, pars), pars(:), sigma, hoptions);
  
  exitflag=0;
  if     strmatch(exitflag, 'tolx')
    exitflag=-5;
    message = [ 'Termination parameter tolerance criteria reached (options.TolX=' ...
              num2str(options.TolX) ')' ];
  elseif strmatch(exitflag, 'tolfun')
    exitflag=-1;
    message = [ 'Termination function tolerance criteria reached (options.TolFun=' ...
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
  constraints = constraints_minmax(pars, constraints);
  [pars,fval,exitflag,output] = GA(@(pars) objective(fun, pars), pars(:),options,constraints);
case {'gradrand','ossrs','fmingradrand'}
% random gradient --------------------------------------------------------------
  [pars,fval,iterations] = ossrs(pars, @(pars) objective(fun, pars), options);
case {'hooke','fminhooke'}
% Hooke-Jeeves direct search ---------------------------------------------------
  [pars,histout] = hooke(pars(:), @(pars) objective(fun, pars), ...
                       options.MaxFunEvals, 2.^(-(0:options.MaxIter)), options.TolFun);
  fval = constraints.criteriaBest;
  iterations      = size(histout,1);
case {'imfil','imfil1','fminimfil'}
% Unconstrained Implicit filtering (version 1998) ------------------------------
  [pars,fval,exitflag,output] = imfil1(pars(:), @(pars) objective(fun, pars), options); 
case {'fminkalman','kalmann','ukfopt'}
% unscented Kalman filter ------------------------------------------------------
  pars = ukfopt(@(pars) objective(fun, pars(:)), pars(:), ...
              options.TolFun, norm(pars)*eye(length(pars)), 1e-6*eye(length(pars)), 1e-6);
  fval = constraints.criteriaBest;
case {'fminlm','LMFsolve'}
% Levenberg-Maquardt steepest descent ------------------------------------------
  if strcmp(options.Display,'iter'), Display = 1;
  else Display = 0; end
  % LMFsolve minimizes the sum of the squares of the objective: sum(objective.^2)
  [pars, fval, iterations, fcount] = LMFsolve(@(pars) objective(fun, pars), pars(:), ...
           'Display',Display, 'FunTol', options.TolFun, 'XTol', options.TolX, ...
           'MaxIter', options.MaxIter, 'Evals',options.MaxFunEvals);
  if iterations < 0
    exitflag=-2;
    iterations = -iterations;
    message='Maximum number of iterations reached';
  elseif fcount < 0
    exitflag=-3;
    message='Maximum number of function evaluations reached'
  else
    exitflag=0;
    message='Algorithm terminated normally'
  end
case {'ntrust','fminnewton'}
% Dogleg trust region, Newton model --------------------------------------------
  [pars,histout,costdata] = ntrust(pars(:),@(pars) objective(fun, pars), ...
       options.TolFun,options.MaxIter);
  fval = constraints.criteriaBest;
  iterations      = size(histout,1);
case {'powell','fminpowell'}
% Powell minimization ----------------------------------------------------------
  if isempty(options.Hybrid), options.Hybrid='Coggins'; end
  if strcmp(lower(options.Hybrid), 'coggins') 
    t = 'Coggins'; method=0;
  else 
    t = 'Golden rule'; method=1;
  end
  options.algorithm  = [ 'Powell Search (by Secchi) [' mfilename '/' t ']' ];
  constraints = constraints_minmax(pars, constraints);
  [pars,fval,exitflag,output] = powell(@(pars) objective(fun, pars), pars, options);
case {'pso','fminpso'}
% particle swarm ---------------------------------------------------------------
  constraints = constraints_minmax(pars, constraints);
  [pars,fval,exitflag,output] = PSO(@(pars) objective(fun, pars),pars(:), ...
     constraints.min(:),constraints.max(:),options);
  message = output.message;
case {'ralg','fminralg','solvopt'}
% Shor's r-algorithm -----------------------------------------------------------
  opt(1) = -1;
  opt(2) = options.TolX;
  opt(3) = options.TolFun;
  opt(4) = options.MaxIter;
  if     strcmp(options.Display,'off') | isempty(options.Display),  opt(5) = -1;
  elseif strcmp(options.Display,'iter'), opt(5) =  1;
  else opt(5)=0; end
  opt(6) = 1e-8; opt(7)=2.5; opt(8)=1e-11;

  % call the optimizer
  [pars,fval,out,iterations, message] = ralg(pars, @(pars) objective(fun, pars), ...
    [], opt,[],[], [], options.MaxFunEvals, options.FunValCheck);
  if out(9) < 0, exitflag = out(9); 
  else exitflag=0; end
case {'fminsearch','fminsearchbnd'}
% Nelder-Mead simplex, with constraints ----------------------------------------
  [pars,fval,exitflag,output] = fminsearch(@(pars) objective(fun, pars), pars, options);
case {'simpsa','fminsimpsa','SIMPSA'}
% simplex/simulated annealing --------------------------------------------------
  constraints = constraints_minmax(pars, constraints);
  [pars,fval,exitflag,output] = SIMPSA(@(pars) objective(fun, pars), pars(:), ...
    constraints.min(:),constraints.max(:),options);
case {'SCE','fminsce'}
% shuffled complex evolution ---------------------------------------------------
  [pars,fval,exitflag,output] = SCE(@(pars) objective(fun, pars), pars(:), ...
    constraints.min(:),constraints.max(:),options);
case {'hPSO','fminswarmhybrid'}
  constraints = constraints_minmax(pars, constraints);
  if isa(options.Hybrid, 'function_handle') | exist(options.Hybrid) == 2
    hoptions.algorithm = [ 'hybrid Particule Swarm Optimizer (by Leontitsis) [' mfilename '/' localChar(options.Hybrid) ']' ];
  else
    hoptions.algorithm = [ 'Particule Swarm Optimizer (by Leontitsis) [fminswarm]' ];
  end

  % transfer optimset options and constraints
  hoptions.space     = [ constraints.min(:) constraints.max(:) ];
  hoptions.MaxIter   = options.MaxIter;
  hoptions.TolFun    = options.TolFun;
  hoptions.TolX      = options.TolX;
  hoptions.Display   = options.Display;
  hoptions.MaxFunEvals=options.MaxFunEvals;
  hoptions.FunValCheck=options.FunValCheck;
  hoptions.OutputFcn  =options.OutputFcn;

  hoptions.Hybrid    = options.Hybrid;
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
    hoptions.maxv = abs(constraints.max(:)-constraints.min(:))/2;
  end
  [pars,fval,exitflag,output] = hPSO(@(pars) objective(fun, pars), pars, hoptions);
otherwise
  options = feval(mfilename, 'defaults');
  [pars,fval,exitflag,output] = fmin_private_wrapper(options.optimizer, fun, pars, ...
    options, constraints, ub);
  return
end

catch
  fval = constraints.criteriaBest;
  pars = constraints.parsBest;
  message = constraints.message;
end

% post optimization checks =====================================================

%pars = apply_constraints(pars, constraints);
if iterations, 
  output.iterations    = iterations;
elseif ~isfield(output,'iterations')
  output.iterations    = constraints.funcCounts;
end
if ~isfield(output,'message')
  output.message         = message;
end
output.funcCount       = constraints.funcCounts;
output.algorithm       = options.algorithm;
output.parsHistory     = constraints.parsHistory;
output.criteriaHistory = constraints.criteriaHistory;
output.parsBest        = constraints.parsBest;
output.criteriaBest    = constraints.criteriaBest;
output.options         = options; 
output.constraints     = constraints;
output.optimizer       = options.optimizer;
output.duration        = etime(clock, t0);

if (exitflag & strcmp(options.Display,'notify')) | ...
   strcmp(options.Display,'final') | strcmp(options.Display,'iter')
  fmin_private_disp_final(output.algorithm, output.message, output.iterations, ...
    output.funcCount, fun, pars, fval);
end

return  % actual end of optimization

% ==============================================================================
% Use nested functions as the criteria wrapper, to access 'constraints' and 'options'
  
  function [c, gc] = objective(fun, pars)
  % criteria to minimize, with gradient support (approx)
  
    % apply constraints on pars first
    pars                = apply_constraints(pars,constraints,options); % private function
    
    % compute criteria
    c = feval(fun, pars);
    if nargout > 1
      gc = gradest(fun, pars(:));
      gc = gc(:);
    end
    
    % check for usual stop conditions MaxFunEvals, TolX, TolFun ..., and call OutputFcn
    [exitflag, message] = fmin_private_std_check(pars, c, ...
       constraints.funcCounts, constraints.funcCounts, options, ...
       constraints.parsPrevious, constraints.criteriaPrevious);
    if exitflag
      constraints.message = message;
      error(message); % will end optimization in try/catch
    end
    
    % save current optimization state
    if c < constraints.criteriaBest, 
      constraints.criteriaBest=c;
      constraints.parsBest    =pars;
    end
    constraints.criteriaPrevious= c;
    constraints.criteriaHistory = [ constraints.criteriaHistory ; c ];
    constraints.parsPrevious    = pars;
    constraints.parsHistory     = [ constraints.parsHistory ; pars ]; 
    constraints.funcCounts      = constraints.funcCounts+1; 
  end

end % optimizer core end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function constraints = constraints_minmax(pars, constraints)
% define default min max in constraints, needed by bounded optimizers
  if ~isfield(constraints, 'min')
    constraints.min = -2*abs(pars); % default min values
    index=find(pars == 0);
    constraints.min(index) = -1;
  end
  if ~isfield(constraints, 'max')
    constraints.max =  2*abs(pars); % default max values
    index=find(pars == 0);
    constraints.max(index) = 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pars,exitflag,message] = apply_constraints(pars, constraints,options)
  % take into account constraints on parameters, and perform stop condition checks
  exitflag=0;
  message='';

  if isfield(constraints, 'fixed')  % fix some parameters
    index = find(constraints.fixed);
    pars(index) = constraints.parsStart(index);
  else
    if isfield(constraints, 'min')  % lower bound for parameters
      index = find(pars(:) < constraints.min(:));
      pars(index) = constraints.min(index);
    end
    if isfield(constraints, 'max')  % upper bound for parameters
      index = find(pars(:) > constraints.max(:));
      pars(index) = constraints.max(index);
    end
    if isfield(constraints, 'step') % restrict parameter change
      parsStep    = pars(:) - constraints.parsPrevious(:);
      index       = find(constraints.steps(:) & abs(parsStep) > abs(constraints.steps(:)) );
      parsStep    = sign(parsStep).*abs(constraints.steps(:));
      pars(index) = constraints.parsPrevious(index) + parsStep(index);
    end
  end
  pars=pars(:)';
end

