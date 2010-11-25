function [pars,fval,exitflag,output] = fmingradrand(varargin)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = FMINGRADRAND(FUN,PARS,[OPTIONS],[CONSTRAINTS]) random gradient optimizer
%
% This minimization method uses a gradient method with random directions.
% Namely, it first determine a random direction in the optimization space
% and then uses a Newton method. This is repeated iteratively until success.
% This method is both fast and less sensitive to local minima traps.
% 
% Calling:
%   fmingradrand(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fmingradrand(fun, pars, options) same as above, with customized options (optimset)
%   fmingradrand(fun, pars, options, fixed) 
%     is used to fix some of the parameters. The 'fixed' vector is then 0 for
%     free parameters, and 1 otherwise.
%   fmingradrand(fun, pars, options, lb, ub) 
%     is used to set the minimal and maximal parameter bounds, as vectors.
%   fmingradrand(fun, pars, options, constraints) 
%     where constraints is a structure (see below).
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fmingradrand(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%     o=fmingradrand('defaults')
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
% Reference: Computer Methods in Applied Mechanics & Engg, Vol  19, (1979) 99
% Contrib: Sheela V. Belur(sbelur@csc.com) 1998 [ossrs]
%
% Version: $Revision: 1.16 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 0 || (nargin == 1 && strcmp(varargin{1},'defaults'))
  options=optimset; % empty structure
  options.Display='';
  options.TolFun =1e-3;
  options.TolX   =1e-8;
  options.MaxIter=1000;
  options.MaxFunEvals=1000;
  options.algorithm  = [ 'Random Gradient (by Belur) [' mfilename ']' ];
  options.optimizer = mfilename;
  pars = options;
  return
end

[pars,fval,exitflag,output] = fmin_private_wrapper(mfilename, varargin{:});

