function [pars,fval,exitflag,output] = fminimfil(varargin)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fminimfil(FUN,PARS,[OPTIONS],[CONSTRAINTS]) Unconstrained Implicit filtering 
%
% Implicit filtering solves unconstrained optimization problems
% Minimization of noisy functions. This is the version 1 of the imfil method.
%
% Calling:
%   fminimfil(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminimfil(fun, pars, options) same as above, with customized options (optimset)
%   fminimfil(fun, pars, options, fixed) 
%     is used to fix some of the parameters. The 'fixed' vector is then 0 for
%     free parameters, and 1 otherwise.
%   fminimfil(fun, pars, options, lb, ub) 
%     is used to set the minimal and maximal parameter bounds, as vectors.
%   fminimfil(fun, pars, options, constraints) 
%     where constraints is a structure (see below).
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminimfil(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%     o=fminimfil('defaults')
%  options.Hybrid specifies the algorithm to use for local hybrid optimizations.
%   This is a string with possible values 'sr1','bfgs' (default),'none'.
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
% Reference: C. T. Kelley, Iterative Methods for Optimization, no. 18 in 
%   Frontiers in Applied Mathematics, SIAM, Philadelphia, 1999.
% Contrib: C. T. Kelley, 1998, Iterative Methods for Optimization [imfil, 1998 version]
%
% Version: $Revision: 1.13 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 0 || (nargin == 1 && strcmp(varargin{1},'defaults'))
  options=optimset;
  % add Matlab std options.
  options.Display='';
  options.TolFun =1e-3;
  options.TolX   =1e-8;
  options.MaxIter=1000;
  options.MaxFunEvals=1000;
  options.Hybrid = 'BFGS';
  options.algorithm  = [ 'Unconstrained Implicit filtering, version 1 (by Kelley) [' mfilename ']' ];
  options.optimizer = mfilename;
  pars = options;
  return
end

[pars,fval,exitflag,output] = fmin_private_wrapper(mfilename, varargin{:});

