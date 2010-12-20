function [pars,fval,exitflag,output] = fmin_private_wrapper(optimizer, fun, pars, options, constraints, ub)
% NOT GOOD: anneal, lm, simpsa
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
% Version: $Revision: 1.5 $
% See also: fminsearch, optimset

% NOTE: all optimizers have been gathered here so that maintenance is minimized
% each user call function only defines the options... The optimozer by itself is
% in the 'private'.
%
% private: 'objective', 'apply_constraints', 'constraints_minmax', 'localChar', 'fmin_private_disp'

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

if strcmp(options.Display,'iter')
  disp([ '** Starting minimization of ' localChar(fun) ' using algorithm ' localChar(options.algorithm) ]);
  spars=pars(1:min(20,length(pars)));
  spars=mat2str(spars');  % as a row
  if length(spars) > 160, spars=[ spars(1:156) ' ...' ]; end
  disp(sprintf('         initial parameters: %s', spars));
  disp('Func_count  min[f(x)]    Parameters');
end

message    = constraints.message;
exitflag   = 0;
iterations = 0;
fval       = Inf;       % in case this is a vector, it should be a row
pars       = pars(:);   % should be a column
output     = [];  
% Optimizer call ===============================================================

try

% calls the optimizer with a wrapped 'objective' function
%    which applies constraints and makes stop condition checks.
% the main optimizer call is within a try/catch block which exits when an early 
%  stop is met. See private 'objective' and 'apply_constraints' below.

switch options.optimizer
case {'fminbfgs','bfgs'}      
% Broyden-Fletcher-Goldfarb-Shanno ---------------------------------------------
  [pars, histout, costdata,iterations] = bfgswopt(pars, @(pars) objective(fun, pars), options.TolFun, options.MaxIter, [], options.TolX);
  iterations = size(histout,1);
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

  [pars, fval, iterations, exitflag, output] = cmaes(@(pars) objective(fun, pars), pars, ...
    sigma, hoptions);

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
  [pars,fval,iretations,output] = GA(@(pars) objective(fun, pars), pars, options,constraints);
case {'gradrand','ossrs','fmingradrand'}
% random gradient --------------------------------------------------------------
  [pars,fval,iterations] = ossrs(pars, @(pars) objective(fun, pars), options);
case {'hooke','fminhooke'}
% Hooke-Jeeves direct search ---------------------------------------------------
  [pars,histout] = hooke(pars, @(pars) objective(fun, pars), ...
                       options.MaxFunEvals, 2.^(-(0:options.MaxIter)), options.TolFun);
  iterations      = size(histout,1);
case {'imfil','fminimfil','fminimfil2','fminimfil3'}
% Unconstrained Implicit filtering (version 1998) ------------------------------
  [pars,fval,iterations,output] = imfil(pars, @(pars) objective(fun, pars), options); 
case {'fminlm','LMFsolve'}
% Levenberg-Maquardt steepest descent ------------------------------------------
  % LMFsolve minimizes the sum of the squares of the objective: sum(objective.^2)
  [pars, fval, iterations, exitflag] = LMFsolve(@(pars) objective_lm(fun, pars), pars, ...
           'Display',0, 'FunTol', options.TolFun, 'XTol', options.TolX, ...
           'MaxIter', options.MaxIter, 'Evals',options.MaxFunEvals);
  switch exitflag
  case -5, message='Termination function tolerance criteria reached';
  case -2, message='Maximum number of iterations reached';
  case -3, message='Maximum number of function evaluations reached';
  case -1, message='Termination parameter tolerance criteria reached';
  end
case {'powell','fminpowell','fminpowell2'}
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
  [pars,fval,exitflag,output] = PSO(@(pars) objective(fun, pars),pars, ...
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
  [pars,fval,exitflag,output] = SIMPSA(@(pars) objective(fun, pars), pars, ...
    constraints.min(:),constraints.max(:),options);
case {'SCE','fminsce'}
% shuffled complex evolution ---------------------------------------------------
  constraints = constraints_minmax(pars, constraints);
  [pars,fval,exitflag,output] = SCE(@(pars) objective(fun, pars), pars, ...
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
  [pars,fval,iterations,output] = hPSO(@(pars) objective(fun, pars), pars, hoptions);
case {'Simplex','fminsimplex'}
% Nelder-Mead simplex state machine --------------------------------------------
  constraints = constraints_minmax(pars, constraints);
  [pars, out]=Simplex('init', pars, abs(constraints.max(:)-constraints.min(:))/10);  % Initialization
  for iterations=1:options.MaxIter
    fval = feval(@(pars) objective(fun, pars), pars);
    [pars,out]=Simplex( fval );
    if Simplex('converged', options.TolFun)             % Test for convergence
      exitflag=-1;
      message= [ 'Termination function tolerance criteria reached (options.TolFun=' ...
                num2str(options.TolFun) ')' ];
      break
    end
    if iterations == options.MaxIter
      exitflag=-1;
      message = [ 'Maximum number of iterations reached (options.MaxIter=' ...
                num2str(options.MaxIter) ')' ];;
      break
    end
  end
  pars=Simplex('centroid'); % obtain the final value.
case {'cgtrust','fmincgtrust'}
% Steihaug Newton-CG-Trust region algoirithm -----------------------------------
  [pars,histout] = cgtrust(pars, @(pars) objective(fun, pars), ...
    [ options.TolFun .1 options.MaxIter options.MaxIter], options.TolX*options.TolX);
% [pars,histout] = levmar(pars, @(pars) objective(fun, pars), options.TolFun, options.MaxIter);
  iterations      = size(histout,1);
otherwise
  options = feval(optimizer, 'defaults');
  [pars,fval,exitflag,output] = fmin_private_wrapper(options.optimizer, fun, pars, ...
    options, constraints, ub);
  return
end % switch
  output.lasterror = lasterror;
end % try


% post optimization checks =====================================================

fval = constraints.criteriaBest;
pars = constraints.parsBest;

if exitflag==0;
  message='Algorithm terminated normally';
end

if iterations, 
  output.iterations    = iterations;
elseif ~isfield(output,'iterations')
  output.iterations    = constraints.funcCounts;
end
if ~isfield(output,'message')
  if isempty(message), message = constraints.message; end
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
  disp([ '** Finishing minimization of ' localChar(fun) ' using algorithm ' localChar(options.algorithm) ]);
  disp(' Func_count     min[f(x)]        Parameters');
  fmin_private_disp(struct('Display','iter'), constraints.funcCounts, fun, pars, fval)
  disp( [ ' Status: ' message ]);
end

return  % actual end of optimization

% ==============================================================================
% Use nested functions as the criteria wrapper, to access 'constraints' and 'options'
  
  function c = objective(fun, pars)
  % criteria to minimize, fun returns a scalar, or vector which is summed
  
    % apply constraints on pars first
    pars                = apply_constraints(pars,constraints,options); % private function
    
    % compute criteria
    c = feval(fun, pars);         % function=row vector, pars=column
    if size(c,1) > 1, c=c'; end
    c  = sum(c(:));
    
    
    % check for usual stop conditions MaxFunEvals, TolX, TolFun ..., and call OutputFcn
    [exitflag, message] = fmin_private_check(pars, c, ...
       constraints.funcCounts, options, constraints.parsPrevious);
    constraints.message = message;
    if exitflag
      error(message); % will end optimization in try/catch
    end
    
    % save current optimization state
    if c < sum(constraints.criteriaBest(:)), 
      constraints.criteriaBest=c(:)';
      constraints.parsBest    =pars;
    end
    constraints.criteriaPrevious= c;
    constraints.criteriaHistory = [ constraints.criteriaHistory ; constraints.criteriaPrevious ];
    constraints.parsPrevious    = pars;
    constraints.parsHistory     = [ constraints.parsHistory ; pars(:)' ]; 
    constraints.funcCounts      = constraints.funcCounts+1; 
    
    if exitflag
      error(message); % will end optimization in try/catch
    end
  end
  
  % LMFsolve supports criteria as a vector of residuals, which sum is the criteria
  % but gradient is used to guide the optimizer
  function c = objective_lm(fun, pars)
  % criteria to minimize, with gradient support (approx)
  
    % apply constraints on pars first
    pars                = apply_constraints(pars,constraints,options); % private function
    
    pars=pars;
    % compute criteria
    c = feval(fun, pars);
    if length(c) == 1,
      c = c*ones(1,10)/10;
    end
    if size(c,1) > 1, c=c'; end
    
    % check for usual stop conditions MaxFunEvals, TolX, TolFun ..., and call OutputFcn
    [exitflag, message] = fmin_private_check(pars, sum(abs(c)), ...
       constraints.funcCounts, options, constraints.parsPrevious);
    constraints.message = message;
    if exitflag
      error(message); % will end optimization in try/catch
    end
    
    % save current optimization state
    if sum(c(:)) < sum(constraints.criteriaBest(:)), 
      constraints.criteriaBest=c(:)';
      constraints.parsBest    =pars;
    end
    constraints.criteriaPrevious= c(:)';
    constraints.criteriaHistory = [ constraints.criteriaHistory ; constraints.criteriaPrevious ];
    constraints.parsPrevious    = pars;
    constraints.parsHistory     = [ constraints.parsHistory ; pars(:)' ]; 
    constraints.funcCounts      = constraints.funcCounts+1; 
    
    if exitflag
      error(message); % will end optimization in try/catch
    end
  end

end % optimizer core end

%   FINJAC       numerical approximation to Jacobi matrix
%   %%%%%%
function J = finjac(FUN,r,x,epsx)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% pars=column, function=row vector or scalar
  lx=length(x);
  J=zeros(lx,length(r));
  if size(x,2) > 1, x=x'; end % column
  if size(r,1) > 1, r=r'; end % row
  if length(epsx)<lx, epsx=epsx*ones(lx,1); end
  for k=1:lx
      dx=.25*epsx(k);
      xd=x;
      xd(k)=xd(k)+dx;
      rd=feval(FUN,xd);
      if size(rd,1) > 1, rd=rd'; end % row
  %   ~~~~~~~~~~~~~~~~    
      if dx, J(k,:)=((rd-r)/dx); end
  end
end

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
  pars=pars';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function strfcn = localChar(fcn)
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

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fmin_private_disp(options, funccount, fun, pars, fval)
% function called during minimization procedure
%
% Displays iteration information every 5 steps, then 10 steps, then 100 steps
% or at every step if iteration is negative


  if funccount > 5
    if     funccount > 500 & mod(funccount,1000) return;
    elseif funccount > 50  & mod(funccount,100) return;
    elseif mod(funccount,10) return; end
  end

  if isfield(options,'Display')
    if strcmp(options.Display, 'iter')
      spars=pars(1:min(20,length(pars)));
      spars=mat2str(spars', 4);  % as a row
      if length(spars) > 50, spars=[ spars(1:47) ' ...' ]; end
      disp(sprintf(' %5.0f    %12.6g   %s', funccount, fval, spars));
    end
  end
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [istop, message] = fmin_private_check(pars, fval, funccount, options, pars_prev, fval_prev)
% standard checks
% fmin_private_check(pars, fval, funccount, options
% or
% fmin_private_check(pars, fval, funccount, options, pars_prev)
% or
% fmin_private_check(pars, fval, funccount, options, pars_prev, fval_prev)
% or
% options=fmin_private_check(options);
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

  % normal terminations: function tolerance reached
  if ~isempty(options.TolFun) && options.TolFun
    if (all(0 < fval) && all(fval <= options.TolFun)) && nargin < 7 % stop on lower threshold
      istop=-1;
      message = [ 'Termination function tolerance criteria reached (fval <= options.TolFun=' ...
                num2str(options.TolFun) ')' ];
    end
    if ~istop & nargin >= 7
      if  all(abs(fval-fval_prev) < options.TolFun*abs(fval)/10) ...
       && all(abs(fval-fval_prev) > 0)
        istop=-12;
        message = [ 'Termination function change tolerance criteria reached (delta(fval) < options.TolFun=' ...
                num2str(options.TolFun) ')' ];
      end
    end
  end
  
  % normal terminations: parameter variation tolerance reached, when function termination is also true
  if (istop==-1 || istop==-12) & nargin >= 6
    if ~isempty(options.TolFun) & options.TolX > 0 ...
      & all(abs(pars-pars_prev(:)) < abs(options.TolX*pars)) ...
      & any(abs(pars-pars_prev(:)) > 0)
      istop=-5;
      message = [ 'Termination parameter tolerance criteria reached (delta(parameters)/parameters <= options.TolX=' ...
            num2str(options.TolX) ')' ];
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
      optimValues.iteration  = iterations;
      optimValues.funcount   = funccount;
      optimValues.fval       = sum(fval(:));
      if isfield(options,'procedure'),        optimValues.procedure=options.procedure;
      elseif isfield(options, 'algorithm'),   optimValues.procedure=options.algorithm;
      else optimValues.procedure  = 'iteration'; end
      istop2 = feval(options.OutputFcn, pars, optimValues, optimValues.state);
      if istop2, 
        istop=-6;
        message = 'Algorithm was terminated by the output function (options.OutputFcn)';
      end
    end
  end
  
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
% -11               Other termination status (cmaes)
% -12               Termination function change tolerance criteria reached

end
