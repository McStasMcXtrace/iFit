function [pars,fval,exitflag,output] = fminmarkov(varargin)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fminmarkov(FUN,PARS,[OPTIONS],[CONSTRAINTS], ...) Markov Chain Monte Carlo optimizer
%
% This minimization method uses a Markov Chain Monte Carlo, optionally 
% with constraints on function parameters.
%
% This is an implementation of the Goodman and Weare 2010 Affine
% invariant ensemble Markov Chain Monte Carlo (MCMC) sampler. MCMC sampling
% enables bayesian inference. The problem with many traditional MCMC samplers
% is that they can have slow convergence for badly scaled problems, and that
% it is difficult to optimize the random walk for high-dimensional problems.
% This is where the GW-algorithm really excels as it is affine invariant. It
% can achieve much better convergence on badly scaled problems. It is much
% simpler to get to work straight out of the box, and for that reason it
% truly deserves to be called the MCMC hammer.
%
% (This code uses a cascaded variant of the Goodman and Weare algorithm).
% 
% Calling:
%   fminmarkov(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminmarkov(fun, pars, options) same as above, with customized options (optimset)
%   fminmarkov(fun, pars, options, fixed) 
%     is used to fix some of the parameters. The 'fixed' vector is then 0 for
%     free parameters, and 1 otherwise.
%   fminmarkov(fun, pars, options, lb, ub) 
%     is used to set the minimal and maximal parameter bounds, as vectors.
%   fminmarkov(fun, pars, options, constraints) 
%     where constraints is a structure (see below).
%   fminmarkov(problem) where problem is a structure with fields
%     problem.objective:   function to minimize
%     problem.x0:          starting parameter values
%     problem.options:     optimizer options (see below)
%     problem.constraints: optimization constraints
%   fminmarkov(..., args, ...)
%     sends additional arguments to the objective function
%       criteria = FUN(pars, args, ...)
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminmarkov(banana,[-1.2, 1])
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
%     o=fminmarkov('defaults')
%  options.MinFunEvals sets the minimum number of function evaluations to reach
%  options.StepSize:  unit-less stepsize (default=2.5).
%  options.ThinChain: Thin all the chains by only storing every N'th step (default=10)
%  options.BurnIn:    fraction of the chain that should be removed. (default=0)
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
% References:
% Goodman & Weare (2010), Ensemble Samplers With Affine Invariance, Comm. App. Math. Comp. Sci., Vol. 5, No. 1, 65-80
% Foreman-Mackey, Hogg, Lang, Goodman (2013), emcee: The MCMC Hammer, arXiv:1202.3665
%
% Contrib:
% By: Aslak Grinsted 2015 https://github.com/grinsted/gwmcmc
%
% Version: $Date$
% See also: fminsearch, optimset

% default options for optimset
if nargin == 0 || (nargin == 1 && strcmp(varargin{1},'defaults'))
  options=optimset; % empty structure
  options.Display='';
  options.TolFun =0;
  options.TolX   =0;
  options.MaxIter=1000;
  options.MaxFunEvals=10000;
  options.PopulationSize=50;
  options.algorithm  = [ 'Cascaded affine invariant ensemble Markov Chain Monte Carlo sampler (by Aslak Grinsted) [' mfilename ']' ];
  options.optimizer = mfilename;
  pars = options;
  return
end

[pars,fval,exitflag,output] = fmin_private_wrapper(mfilename, varargin{:});

