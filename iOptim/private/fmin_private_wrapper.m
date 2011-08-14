function [pars,fval,exitflag,output] = fmin_private_wrapper(optimizer, fun, pars, options, constraints, ub)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fmin_private_wrapper(OPTIMIZER, FUN,PARS,[OPTIONS],[constraints]) wrapper to optimizers
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
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fmin_private_wrapper('fminimfil',banana,[-1.2, 1])
%
% Input:
%  OPTIMIZER is the name/handle to an optimizer function, or '' for default
%  FUN is a function handle (anonymous function or inline) with a loss
%  function, which may be of any type, and needn't be continuous. It does,
%  however, need to return a single/vector value.
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the simulated annealing, 
%  compliant with optimset. Default options may be obtained with
%     o=fmin_private_wrapper(optimizer,'defaults')
%
%  constraints may be specified as a structure
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
% Version: $Revision: 1.23 $
% See also: fminsearch, optimset

% NOTE: all optimizers have been gathered here so that maintenance is minimized
% each user call function only defines the options... The optimizer by itself is
% in the 'private'.
%
% private: 'inline_objective', 
%          'inline_apply_constraints', 
%          'inline_constraints_minmax', 
%          'inline_localChar', 
%          'inline_disp'
%          'inline_estimate_uncertainty'
%

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

% parameter handling ===========================================================

% nargin stuff (number of parameters)
% default options for optimset

if nargin < 1, optimizer=''; end
if nargin < 2, fun = ''; end
if isempty(optimizer), optimizer = 'fminimfil'; end
if isempty(fun),       fun = 'defaults'; end

if nargin < 4
  options=[];
end
if nargin < 5
  constraints = [];
end
if nargin < 6
  ub = [];
end

if strcmp(fun,'defaults')
  pars = feval(optimizer, 'defaults');
  return
elseif nargin < 2
  error([ 'syntax is: ' optimizer '(optimizer, inline_objective, parameters, ...)' ] );
elseif nargin == 2 && isstruct(fun)
  if     isfield(fun, 'x0'),          pars=fun.x0;
  elseif isfield(fun, 'guess'),       pars=fun.guess;
  elseif isfield(fun, 'Guess'),       pars=fun.Guess; end
  if     isfield(fun, 'options'),     options=fun.options; end
  if     isfield(fun, 'constraints'), constraints=fun.constraints; end 
  if     isfield(fun, 'inline_objective'),   tmp=fun.inline_objective; fun=[]; fun=tmp; 
  elseif isfield(fun, 'model'),       tmp=fun.model; fun=[]; fun=tmp;
  elseif isfield(fun, 'f'),           tmp=fun.f; fun=[]; fun=tmp;
  elseif isfield(fun, 'function'),    tmp=fun.function; fun=[]; fun=tmp; end
elseif nargin < 3
  error([ 'syntax is: ' inline_localChar(optimizer) '(inline_objective, parameters, ...)' ] );
end

% default arguments when missing
if isempty(options)
  options=feval(optimizer, 'defaults');
end
if length(constraints) & isnumeric(constraints) 
  if nargin == 4,               % given as fixed index vector
    fixed = constraints; constraints=[];
    constraints.fixed = fixed;  % avoid warning for variable redefinition.
  elseif nargin == 5            % given as lb,ub parameters (nargin==5)
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
constraints.funcCount       = 0;
constraints.message         = '';

options.optimizer = optimizer;

options=fmin_private_check(options, feval(options.optimizer,'defaults'));
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

if ischar(options.TolFun)
  options.TolFunChar = options.TolFun;
  if options.TolFun(end)=='%'
    options.TolFun(end)='';
    fval = inline_objective(fun, pars);
    options.TolFun = abs(str2num(options.TolFun)*fval/100);
  else
    options.TolFun = str2num(options.TolFun);
  end
else
  fval = NaN;
end

if ischar(options.TolX)
  options.TolXChar=options.TolX;
  options.TolX = abs(str2num(options.TolX)*pars(:)/100);
end

if strcmp(options.Display,'iter')
  disp([ '** Starting minimization of ' inline_localChar(fun) ' using algorithm ' inline_localChar(options.algorithm) ]);
  disp('Func_count  min[f(x)]    Parameters');
  inline_disp(options, constraints.funcCount , fun, pars, fval)
end

message    = constraints.message;
exitflag   = 0;
iterations = 0;
fval       = Inf;       % in case this is a vector, it should be a row
pars       = pars(:);   % should be a column
output     = [];  
% Optimizer call ===============================================================

try

% calls the optimizer with a wrapped 'inline_objective' function
%    which applies constraints and makes stop condition checks.
% the main optimizer call is within a try/catch block which exits when an early 
%  stop is met. See private 'inline_objective' and 'inline_apply_constraints' below.

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

  [pars, fval, iterations, exitflag, output] = cmaes(@(pars) inline_objective(fun, pars), pars, ...
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
  constraints = inline_constraints_minmax(pars, constraints);
  [pars,fval,iretations,output] = GA(@(pars) inline_objective(fun, pars), pars, options,constraints);
case {'gradrand','ossrs','fmingradrand'}
% random gradient --------------------------------------------------------------
  [pars,fval,iterations] = ossrs(pars, @(pars) inline_objective(fun, pars), options);
case {'hooke','fminhooke'}
% Hooke-Jeeves direct search ---------------------------------------------------
  [pars,histout] = hooke(pars, @(pars) inline_objective(fun, pars), ...
                       options.MaxFunEvals, 2.^(-(0:options.MaxIter)), options.TolFun);
  iterations      = size(histout,1);
case {'imfil','fminimfil'}
% Unconstrained Implicit filtering (version 1998) ------------------------------
  [pars,fval,iterations,output] = imfil(pars, @(pars) inline_objective(fun, pars), options); 
case {'fminlm','LMFsolve'}
% Levenberg-Maquardt steepest descent ------------------------------------------
  % LMFsolve minimizes the sum of the squares of the inline_objective: sum(inline_objective.^2)
  [pars, fval, iterations, exitflag] = LMFsolve(@(pars) inline_objective_lm(fun, pars), pars, ...
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
    t = 'Coggins'; method=0;
  else 
    t = 'Golden rule'; method=1;
  end
  options.algorithm  = [ 'Powell Search (by Secchi) [' options.optimizer '/' t ']' ];
  constraints = inline_constraints_minmax(pars, constraints);
  [pars,fval,exitflag,output] = powell(@(pars) inline_objective(fun, pars), pars, options);
case {'pso','fminpso'}
% particle swarm ---------------------------------------------------------------
  constraints = inline_constraints_minmax(pars, constraints);
  [pars,fval,exitflag,output] = PSO(@(pars) inline_objective(fun, pars),pars, ...
     constraints.min(:),constraints.max(:),options);
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
  [pars,fval,out,iterations, message] = ralg(pars, @(pars) inline_objective(fun, pars), ...
    [], opt,[],[], [], options.MaxFunEvals, options.FunValCheck);
  if out(9) < 0, exitflag = out(9); 
  else exitflag=0; end
case {'fminsearch','fminsearchbnd'}
% Nelder-Mead simplex, with constraints ----------------------------------------
  [pars,fval,exitflag,output] = fminsearch(@(pars) inline_objective(fun, pars), pars, options);
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
  constraints = inline_constraints_minmax(pars, constraints);
  [pars,fval,exitflag,output] = SIMPSA(@(pars) inline_objective(fun, pars), pars, ...
    constraints.min(:),constraints.max(:),options);
case {'SCE','fminsce'}
% shuffled complex evolution ---------------------------------------------------
  constraints = inline_constraints_minmax(pars, constraints);
  [pars,fval,exitflag,output] = SCE(@(pars) inline_objective(fun, pars), pars, ...
    constraints.min(:),constraints.max(:),options);
case {'hPSO','fminswarmhybrid'}
  constraints = inline_constraints_minmax(pars, constraints);
  if isa(options.Hybrid, 'function_handle') | exist(options.Hybrid) == 2
    hoptions.algorithm = [ 'hybrid Particle Swarm Optimizer (by Leontitsis) [' options.optimizer '/' inline_localChar(options.Hybrid) ']' ];
  else
    hoptions.algorithm = [ 'Particle Swarm Optimizer (by Leontitsis) [fminswarm]' ];
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
  [pars,fval,iterations,output] = hPSO(@(pars) inline_objective(fun, pars), pars, hoptions);
case {'Simplex','fminsimplex'}
% Nelder-Mead simplex state machine --------------------------------------------
  constraints = inline_constraints_minmax(pars, constraints);
  [pars, out]=Simplex('init', pars, abs(constraints.max(:)-constraints.min(:))/10);  % Initialization
  for iterations=1:options.MaxIter
    fval = feval(@(pars) inline_objective(fun, pars), pars);
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
  [pars,histout] = cgtrust(pars, @(pars) inline_objective(fun, pars), ...
    [ options.TolFun .1 options.MaxIter options.MaxIter], options.TolX*options.TolX);
% [pars,histout] = levmar(pars, @(pars) inline_objective(fun, pars), options.TolFun, options.MaxIter);
  iterations      = size(histout,1);
% not so efficient optimizers ==================================================
case {'fminanneal','anneal'}  
% simulated annealing ----------------------------------------------------------
  options.MaxTries   = options.MaxIter/10;
  options.StopVal    = options.TolFun;
  options.Verbosity=0;
  [pars,fval,iterations,exitflag] = anneal(@(pars) inline_objective(fun, pars), pars(:)', options);
  if exitflag==-7,message='Maximum consecutive rejections exceeded (anneal)'; end
case {'fminbfgs','bfgs'}      
% Broyden-Fletcher-Goldfarb-Shanno ---------------------------------------------
  [pars, histout, costdata,iterations] = bfgswopt(pars(:), @(pars) inline_objective(fun, pars), options.TolFun, options.MaxIter);
  iterations = size(histout,1);
case {'fminkalman','kalmann','ukfopt'}
% unscented Kalman filter ------------------------------------------------------
  [pars,iterations] = ukfopt(@(pars) inline_objective(fun, pars(:)), pars(:), ...
              options.TolFun, norm(pars)*eye(length(pars)), 1e-6*eye(length(pars)), 1e-6);
case {'ntrust','fminnewton'}
% Dogleg trust region, Newton model --------------------------------------------
  [pars,histout,costdata] = ntrust(pars(:),@(pars) inline_objective(fun, pars), ...
       options.TolFun,options.MaxIter);
  iterations      = size(histout,1);
case {'buscarnd','fminrand'}
% adaptive random search -------------------------------------------------------
  [pars,fval]=buscarnd(@(pars) inline_objective(fun, pars), pars, options);
otherwise
  options = feval(optimizer, 'defaults');
  [pars,fval,exitflag,output] = fmin_private_wrapper(options.optimizer, fun, pars, ...
    options, constraints, ub);
  return
end % switch
catch
  output.lasterror = lasterror;
end % try

if isstruct(output) && isfield(output,'lasterror') && isempty(strfind(output.lasterror.message, 'stop condition:'))
  disp('Code error when launching the optimizer. Please fix it...')
  disp(output.lasterror.message);
  for index=1:length(output.lasterror.stack)
    disp(output.lasterror.stack(index))
  end
  rethrow(output.lasterror);
end

% post optimization checks =====================================================

fval = constraints.criteriaBest;
pars = constraints.parsBest;

% raise fminplot if is exists
if ~isempty(options.OutputFcn) & strcmp(options.OutputFcn, 'fminplot')
  h = findall(0, 'Tag', 'fminplot');
  if ~isempty(h), 
    figure(h(1)); 
    set(h, 'Visible', 'on');
  end
end

if exitflag==0;
  message='Algorithm terminated';
end

if iterations, 
  output.iterations    = iterations;
elseif ~isfield(output,'iterations')
  output.iterations    = constraints.funcCount ;
end
if ~isfield(output,'message')
  if isempty(message), message = constraints.message; end
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

% estimate parameter uncertainty from the search trajectory
index      = find(output.criteriaHistory < min(output.criteriaHistory)*4);   % identify tolerance region around optimum                       
delta_pars = (output.parsHistory(index,:)-repmat(output.parsBest,[length(index) 1])); % get the corresponding parameter set
weight_pars= exp(-((output.criteriaHistory(index)-min(output.criteriaHistory))/min(output.criteriaHistory)).^2 / 8); % Gaussian weighting for the parameter set
weight_pars= repmat(weight_pars,[1 length(output.parsBest)]);
output.parsHistoryUncertainty = sqrt(sum(delta_pars.*delta_pars.*weight_pars)./sum(weight_pars));

if length(pars)^2*output.fevalDuration/2 < 60 % should spend less than a minute to compute the Hessian
  try
  [dp, covp, corp,jac,hessian]  = inline_estimate_uncertainty(fun, pars);
  if ~isempty(covp)
    output.parsHessianUncertainty = reshape(abs(dp), size(pars));
    output.parsHessianCovariance  = covp;
    output.parsHessianCorrelation = corp;
    output.parsHessian            = hessian;
    output.parsJacobian           = jac;
  end
  end
else
  output.parsHessianUncertainty = [];
  output.parsHessianCovariance  = [];
  output.parsHessianCorrelation = [];
  output.parsHessian            = [];
  output.parsJacobian           = [];
end

if strcmp(options.Display,'final') | strcmp(options.Display,'iter')
  disp([ sprintf('\n') '** Finishing minimization of ' inline_localChar(fun) ' using algorithm ' inline_localChar(options.algorithm) ]);
  disp( [ ' Status: ' message ]);
  disp(' Func_count     min[f(x)]        Parameters');
  inline_disp(struct('Display','iter'), -constraints.funcCount , fun, pars, mean(fval));
  
  if length(index) > 10
    disp(' Gaussian uncertainty on parameters (half width, from the optimization history)')
    inline_disp(struct('Display','iter'), -1, fun, output.parsHistoryUncertainty, NaN);
  else
    disp(' Gaussian uncertainty on parameters (half width, from the Hessian matrix)')
    inline_disp(struct('Display','iter'), -1, fun, output.parsHessianUncertainty, NaN);
  end
end

return  % actual end of optimization

% ==============================================================================
% Use nested functions as the criteria wrapper, to access 'constraints' and 'options'
  
  function c = inline_objective(fun, pars)
  % criteria to minimize, fun returns a scalar, or vector which is summed
  
    % apply constraints on pars first
    pars                = inline_apply_constraints(pars,constraints,options); % private function
    
    % compute criteria
    t = clock;
    c = feval(fun, pars);         % function=row vector, pars=column
    c  = sum(c);
    
    % check for usual stop conditions MaxFunEvals, TolX, TolFun ..., and call OutputFcn
    [exitflag, message] = fmin_private_check(pars, c, ...
       constraints.funcCount , options, constraints);
    constraints.message = message;
    
    % save current optimization state
    if c < sum(constraints.criteriaBest(:)), 
      constraints.criteriaBest=c;
      constraints.parsBest    =pars;
    end
    constraints.fevalDuration   = etime(clock, t); % time required to estimate the criteria
    constraints.criteriaPrevious= c;
    constraints.criteriaHistory = [ constraints.criteriaHistory ; sum(constraints.criteriaPrevious) ];
    constraints.parsPrevious    = pars;
    constraints.parsHistory     = [ constraints.parsHistory ; pars ]; 
    constraints.funcCount       = constraints.funcCount +1; 
    
    if exitflag
      error([ 'stop condition: ' message ]); % will end optimization in try/catch
    end
  end
  
  % LMFsolve supports criteria as a vector of residuals, which sum is the criteria
  % but gradient is used to guide the optimizer
  function c = inline_objective_lm(fun, pars)
  % criteria to minimize, with gradient support (approx)
    
    % apply constraints on pars first
    pars                = inline_apply_constraints(pars,constraints,options); % private function

    % compute criteria
    t = clock;
    c = feval(fun, pars); % function=row vector, pars=column
    if length(c) == 1,
      c = c*ones(1,10)/10;
    end
    if length(c) > 1, c=c(:)'; end
    
    % check for usual stop conditions MaxFunEvals, TolX, TolFun ..., and call OutputFcn
    [exitflag, message] = fmin_private_check(pars, sum(abs(c)), ...
       constraints.funcCount , options, constraints);
    constraints.message = message;
    
    % save current optimization state
    if sum(c) < sum(constraints.criteriaBest), 
      constraints.criteriaBest=c;
      constraints.parsBest    =pars;
    end
    constraints.fevalDuration   = etime(clock, t); % time required to estimate the criteria
    constraints.criteriaPrevious= sum(c);
    constraints.criteriaHistory = [ constraints.criteriaHistory ; constraints.criteriaPrevious ];
    constraints.parsPrevious    = pars;
    constraints.parsHistory     = [ constraints.parsHistory ; pars ]; 
    constraints.funcCount       = constraints.funcCount +1; 
    
    if exitflag
      error([ 'stop condition: ' message ]); % will end optimization in try/catch
    end
  end

end % fmin_private_wrapper optimizer core end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function constraints = inline_constraints_minmax(pars, constraints)
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
end % inline_constraints_minmax

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pars,exitflag,message] = inline_apply_constraints(pars, constraints,options)
  % take into account constraints on parameters, and perform stop condition checks
  exitflag=0;
  message='';

  if isfield(constraints, 'fixed')  % fix some parameters
    index = find(constraints.fixed & ~isnan(constraints.fixed));
    if ~isempty(index), pars(index) = constraints.parsStart(index); end
  else
    if isfield(constraints, 'min')  % lower bound for parameters
      index = find(pars(:) < constraints.min(:) & ~isnan(constraints.min(:)));
      if ~isempty(index), pars(index) = constraints.min(index); end
    end
    if isfield(constraints, 'max')  % upper bound for parameters
      index = find(pars(:) > constraints.max(:) & ~isnan(constraints.max(:)));
      if ~isempty(index), pars(index) = constraints.max(index); end
    end
    if isfield(constraints, 'step') % restrict parameter change
      parsStep    = pars(:) - constraints.parsPrevious(:);
      index       = find(constraints.steps(:) & abs(parsStep) > abs(constraints.steps(:)) & ~isnan(constraints.steps(:)));
      if ~isempty(index), 
        parsStep    = sign(parsStep).*abs(constraints.steps(:));
        pars(index) = constraints.parsPrevious(index) + parsStep(index);
      end
    end
  end
  pars=pars(:)';
end % inline_apply_constraints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function strfcn = inline_localChar(fcn)
% Convert the fcn to a string for printing

  if ischar(fcn)
      strfcn = fcn;
  elseif isa(fcn,'inline')
      strfcn = char(fcn);
  elseif isa(fcn,'function_handle')
      strfcn = func2str(fcn);
  else
      try
          strfcn = char(fcn);
      catch
          strfcn = '(name not printable)';
      end
  end

end % inline_localChar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inline_disp(options, funccount, fun, pars, fval)
% function called during minimization procedure
%
% Displays iteration information every 5 steps, then 10 steps, then 100 steps
% or at every step if iteration is negative


  if funccount > 5
    if funccount > 50  & mod(funccount,100) return;
    elseif mod(funccount,10) return; end
  end

  if isfield(options,'Display')
    if strcmp(options.Display, 'iter')
      spars=pars(1:min(20,length(pars)));
      spars=mat2str(spars', 4);  % as a row
      if length(spars) > 50, spars=[ spars(1:47) ' ...' ]; end
      if isfinite(funccount) && isfinite(fval)
        disp(sprintf(' %5.0f    %12.6g          %s', abs(funccount), sum(fval), spars));
      else
        disp(sprintf('                                %s', spars));
      end
    end
  end
  
end % inline_disp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [istop, message] = fmin_private_check(pars, fval, funccount, options, constraints)
% standard checks
% fmin_private_check(pars, fval, funccount, options, constraints)
% or
% options=fmin_private_check(options, default_options);

  istop=0; message='';
  
  % check of option members
  if nargin<=2
    options=pars;
    if nargin ==2, 
      default=fval; 
      checks=fieldnames(default);
    else 
      fval=[]; 
      checks={'TolFun','TolX','Display','MaxIter','MaxFunEvals','FunValCheck','OutputFcn','algorithm'};
    end
    
    for index=1:length(checks)
      if ~isfield(options, checks{index}), 
        if isfield(default, checks{index}), 
          options=setfield(options,checks{index},getfield(default, checks{index}));
        else
          options=setfield(options,checks{index},[]); 
        end
      end
    end
    istop=options;
    return
  end
  
  pars_prev = constraints.parsPrevious;
  fval_prev = constraints.criteriaPrevious;
  fval_best = constraints.criteriaBest;
  fval_mean = mean(constraints.criteriaHistory(max(1, length(constraints.criteriaHistory)-10):end));
  
  % handle relative stop conditions
  if isfield(options,'TolFunChar')
    options.TolFun = options.TolFunChar;
  end
  if ischar(options.TolFun)
    if options.TolFun(end)=='%'
      options.TolFun(end)='';
      options.TolFun = abs(str2num(options.TolFun)*fval/100);
    else
      options.TolFun = str2num(options.TolFun);
    end
  end
  if isfield(options,'TolXChar')
    options.TolX = options.TolXChar;
  end
  if ischar(options.TolX)
    if options.TolX(end)=='%'
      options.TolX(end)='';
      options.TolX = abs(str2num(options.TolX)*pars(:)/100);
    else
      options.TolX = str2num(options.TolX);
    end
  end

  % normal terminations: function tolerance reached
  if ~isempty(options.TolFun) && options.TolFun ~= 0 && funccount >= 5*length(pars)
    if (all(0 < fval) && all(fval <= options.TolFun)) % stop on lower threshold
      istop=-1;
      message = [ 'Converged: Termination function tolerance criteria reached (fval <= options.TolFun=' ...
                num2str(options.TolFun) ')' ];
    end
    if ~istop
      % stop on criteria change
      if  all(abs(fval-fval_prev) < options.TolFun) ...
       && all(abs(fval-fval_prev) > 0) ...
       && all(fval < fval_mean - options.TolFun) 
        istop=-12;
        message = [ 'Converged: Termination function change tolerance criteria reached (delta(fval) < options.TolFun=' ...
                num2str(options.TolFun) ')' ];
      end
    end
  end
  
  % normal terminations: parameter variation tolerance reached, when function termination is also true
  if (istop==-1 || istop==-12) 
    if ~isempty(options.TolX) && options.TolX > 0 ...
      && all(abs(pars(:)-pars_prev(:)) < abs(options.TolX)) ...
      && any(abs(pars(:)-pars_prev(:)) > 0)
      istop=-5;
      message = [ 'Converged: Termination parameter tolerance criteria reached (delta(parameters) <= options.TolX=' ...
            num2str(mean(options.TolX)) ')' ];
    end
  end
  
  % abnormal terminations
  if ~istop

    if options.MaxFunEvals > 0 & funccount >= options.MaxFunEvals
      istop=-3;
      message = [ 'Maximum number of function evaluations reached (options.MaxFunEvals=' ...
                num2str(options.MaxFunEvals) ')' ];
    end
    
    if strcmp(options.FunValCheck,'on') & any(isnan(fval) | isinf(fval))
      istop=-4;
      message = 'Function value is Inf or Nan (options.FunValCheck)';
    end

    if ~isempty(options.OutputFcn)
      optimValues = options;
      if ~isfield(optimValues,'state')
        if istop,               optimValues.state='done';
        elseif funccount  <= 5, optimValues.state='init';
        else                    optimValues.state='iter'; end
      end
      optimValues.funcount   = funccount;
      optimValues.funcCount  = funccount;
      optimValues.funccount  = funccount;
      optimValues.fval       = sum(fval(:));
      if isfield(options,'procedure'),        optimValues.procedure=options.procedure;
      elseif isfield(options, 'algorithm'),   optimValues.procedure=options.algorithm;
      elseif isfield(options, 'optimizer'),   optimValues.procedure=options.optimizer;
      else optimValues.procedure  = 'iteration'; end
      istop2 = feval(options.OutputFcn, pars, optimValues, optimValues.state);
      if istop2, 
        istop=-6;
        message = 'Algorithm was terminated by the output function (options.OutputFcn)';
      end
    end
    if istop
      funccount = -funccount; % trigger iteration display
    end
    inline_disp(options,  funccount, options.optimizer, pars, fval);
  end

end % fmin_private_check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dp,covp,corp,jac,hessian] = inline_estimate_uncertainty(fun, pars, options)
% [dp,covp,corp] = inline_estimate_uncertainty(fun, pars, options)
% 
% Estimates the uncertainty around an optimization solution using
% the error matrix from the criteria jacobian inversion.
% 
% Calling:
%   [p,c,e,o]  = fminpso(fun, p0);
%   dp = o.parsHessianUncertainty;
%
% Input:
%  FUN is a function handle (anonymous function or inline) with a loss
%  function, which may be of any type, and needn't be continuous. It does,
%  however, need to return a single value.
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the simulated annealing, 
%  compliant with optimset. Default options may be obtained with
%     o=fmin_private_wrapper(optimizer,'defaults')
%
% Output:
%  DP is the gaussian uncertainty around PARS
%  COVP is the error matrix 
%  CORP is the correlation matrix
%  JAC  is the Jacobian
%  HESSIAN is the Hessian


  n=length(pars);
  if nargin < 3, options=[]; end
  if isfield(options,'TolX') TolX = options.TolX; 
  else TolX = 0; end

  % initialize the curvature matrix alpha = '1/2 d2 Chi2/dpi/dpj' (half Hessian)
  alpha= zeros(n);
  dp   = zeros(size(pars));
  chisq= sum(feval(fun, pars));
  
  covp = [];
  corp = [];
  jac  = [];
  hessian = [];
  if TolX <= 0, 
    TolX = 0.01*pars;
  end
  if length(TolX) == 1
    dp   = TolX*ones(size(pars));
  else
    dp   = TolX;
  end

  % we now build the error matrix 'alpha' and the Jacobian
  jac = zeros(n,1);
  for i=1:n
    p    = pars; p(i) = p(i)+dp(i); chi1 = sum(feval(fun, p));
    p    = pars; p(i) = p(i)-dp(i); chi2 = sum(feval(fun, p));
    alpha(i,i) = (chi1-2*chisq+chi2)/2/dp(i)/dp(i); % diagonal terms
    jac(i) = (chi1-chisq)/dp(i);
    
    for j=i+1:n
      p=pars; p(i)=p(i)+dp(i); p(j)=p(j)+dp(j); chi1=sum(feval(fun,p));
      p=pars; p(i)=p(i)+dp(i); p(j)=p(j)-dp(j); chi2=sum(feval(fun,p));
      p=pars; p(i)=p(i)-dp(i); p(j)=p(j)+dp(j); chi3=sum(feval(fun,p));
      p=pars; p(i)=p(i)-dp(i); p(j)=p(j)-dp(j); chi4=sum(feval(fun,p));
      alpha(i,j)=(chi1-chi2-chi3+chi4)/8/dp(i)/dp(j);
      alpha(j,i)=alpha(i,j); % off diagonal terms (symmetric)
    end
  end
  if any(isnan(alpha(:))), return; end 
  hessian=2*alpha;
  alpha = alpha/chisq;      % normalized error matrix
  covp  = pinv(alpha);       % COV MATRIX
  dp    = sqrt(abs(diag(covp))); % uncertainty on parameters
  corp  = covp./(dp*dp');   % correlation matrix
  
end % inline_estimate_uncertainty

