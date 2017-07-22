function [pars,fval,exitflag,output] = fmin(varargin)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = FMIN(FUN,PARS,[OPTIONS],[CONSTRAINTS], ...) Best optimizer
%
% This minimization method is determined automatically from the objective function
% behaviour and number of free parameters. You can however force a specific 
% optimizer by setting e.g. options.optimizer='fminpso'
%
% WARNING: as the selected optimizer may change from one call to an other, the
% solution found may vary as well. To avoid that, rather use a specific optimizer.
%
% Best optimizers are:
%   fminpso:    Particle Swarm Optimization
%   fminpowell: Powell with Coggins line search
%   fminhooke:  Hooke-Jeeves direct search
%   fminralg:   Shor R-algorithm
%   fminsimpsa: Simplex/simulated annealing
%   fminimfil:  Unconstrained Implicit filtering
% Type <a href="matlab:doc(iData,'Optimizers')">doc(iData,'Optimizers')</a> to access the Optimizers Documentation.
%
% Calling:
%   fmin(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fmin(fun, pars, options) same as above, with customized options (optimset)
%   fmin(fun, pars, options, fixed) 
%     is used to fix some of the parameters. The 'fixed' vector is then 0 for
%     free parameters, and 1 otherwise.
%   fmin(fun, pars, options, lb, ub) 
%     is used to set the minimal and maximal parameter bounds, as vectors.
%   fmin(fun, pars, options, constraints) 
%     where constraints is a structure (see below).
%   fmin(problem) where problem is a structure with fields
%     problem.objective:   function to minimize
%     problem.x0:          starting parameter values
%     problem.options:     optimizer options (see below)
%     problem.constraints: optimization constraints
%   fmin(..., args, ...)
%     sends additional arguments to the objective function
%       criteria = FUN(pars, args, ...)
%
% The options structure may contain the following members, in agreement with 'optimset':
%    options.Display: Level of display [ off | iter | notify | final ]. Default is 'off'
%    options.MaxFunEvals: Maximum number of function evaluations allowed, sometimes referred as the 'cost' or 'budget'.
%    options.MaxIter: Maximum number of iterations allowed
%    options.TolFun: Termination tolerance on the function value (absolute value or change). Use 'x%' to specify a relative function change.
%    options.TolX: Termination tolerance on parameter change. Use 'x%' to specify a relative parameter change.
%    options.OutputFcn: Name of an output function. When set, it is called at each iteration step. You may use 'fminplot', which is provided in Optimizers. Refer to the Fit page for more information about fminplot. A simpler/faster alternative is the 'fminstop' option.
%    options.PlotFcns: same as OutputFcn, but can be a set of function in a cell array.
%    options.FunValCheck: Check for invalid values, such as NaN or complex
%    options.MinFunEvals: when set, waits for a given number of iterations before testing for convergence
%    options.optimizer: the optimizer to use
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fmin(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string): criteria = FUN(PARS)
%  It needs to return a single value or vector.
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess. PARS can also be given as a single-level structure.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%     o=fmin('defaults')
%  options.MinFunEvals sets the minimum number of function evaluations to reach
%  An empty OPTIONS sets the default configuration.
%
%  CONSTRAINTS may be specified as a structure
%   constraints.min= vector of minimal values for parameters
%   constraints.max= vector of maximal values for parameters
%   constraints.fixed= vector having 0 where parameters are free, 1 otherwise
%   constraints.step=  vector of maximal parameter changes per iteration
%   constraints.eval=  expression making use of 'p', 'constraints', and 'options' 
%                        and returning modified 'p'
%                      or function handle p=@constraints.eval(p)
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
% See also: fminsearch, optimset, 
% (c) E.Farhi, ILL. License: EUPL.

% default options for optimset
if nargin == 0 || (nargin == 1 && strcmp(varargin{1},'defaults'))
  options=optimset; % empty structure
  options.Display='';
  options.TolFun =1e-3;
  options.TolX   =1e-8;
  options.MaxIter=1000;
  options.MaxFunEvals=10000;
  options.PopulationSize=50;
  options.algorithm  = [ 'Automatic optimizer [' mfilename ']' ];
  options.optimizer = mfilename;
  pars = options;
  return
end

[pars,fval,exitflag,output] = fmin_private_wrapper(mfilename, varargin{:});

