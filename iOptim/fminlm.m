function [pars,fval,exitflag,output] = fminlm(varargin)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fminlm(FUN,PARS,[OPTIONS],[CONSTRAINTS]) Levenberg-Maquardt search
%
% This minimization method uses the Levenberg-Maquardt steepest descent 
% in Least-Squares Sense. It finds parameters in order to bring the objective
% to zero (and not to its lowest value). This implementation is not as efficient
% as when used directly with residuals.
% 
% Calling:
%   fminlm(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminlm(fun, pars, options) same as above, with customized options (optimset)
%   fminlm(fun, pars, options, fixed) 
%     is used to fix some of the parameters. The 'fixed' vector is then 0 for
%     free parameters, and 1 otherwise.
%   fminlm(fun, pars, options, lb, ub) 
%     is used to set the minimal and maximal parameter bounds, as vectors.
%   fminlm(fun, pars, options, constraints) 
%     where constraints is a structure (see below).
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminlm(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%     o=fminlm('defaults')
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
% Reference: 
% Fletcher, R., (1971) Rpt. AERE-R 6799, Harwell
% Fletcher, R., Computer Journal 1970, 13, 317-322
% Contrib: Miroslav Balda, balda AT cdm DOT cas DOT cz 2009 [LMFsolve]
%
% Version: $Revision: 1.9 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 0 || (nargin == 1 && strcmp(varargin{1},'defaults'))
  options=optimset; % empty structure
  options.Display  = [];        %   no print of iterations
  options.MaxIter  = 5000;       %   maximum number of iterations allowed
  options.ScaleD   = [];        %   automatic scaling by D = diag(diag(J'*J))
  options.TolFun   = 1e-5;      %   tolerace for final function value
  options.TolX     = 1e-4;      %   tolerance on difference of x-solutions
  options.MaxFunEvals=10000;
  options.algorithm  = [ 'Levenberg-Maquardt (by Balda) [' mfilename ']' ];
  options.optimizer = mfilename;
  pars = options;
  return
end

[pars,fval,exitflag,output] = fmin_private_wrapper(mfilename, varargin{:});

