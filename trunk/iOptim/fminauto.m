function [pars,fval,exitflag,output] = fminauto(fun, pars, options, varargin)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = FMINAUTO(FUN,PARS,[OPTIONS],[CONSTRAINTS]) Automatic best Optimization
%
% This minimization method uses a set of algorithms for 
% finding the minimum of the function 'FUN' in the real space. 
%
% Calling:
%   fminauto(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminauto(fun, pars, options) same as above, with customized options (optimset)
%   fminauto(fun, pars, options, fixed) 
%     is used to fix some of the parameters. The 'fixed' vector is then 0 for
%     free parameters, and 1 otherwise.
%   fminauto(fun, pars, options, lb, ub) 
%     is used to set the minimal and maximal parameter bounds, as vectors.
%   fminauto(fun, pars, options, constraints) 
%     where constraints is a structure (see below).
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminauto(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%      o=fminauto('defaults');
%
%  CONSTRAINTS may be specified as a structure
%   constraints.min= vector of minimal values for parameters
%   constraints.max= vector of maximal values for parameters
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
% See also: fminsearch, optimset


if nargin == 1 & strcmp(fun,'defaults')
  options=fminswarmhybrid('defaults');
  options.Hybrid='none';
  options.TolX   =1e-12;
  options.algorithm = [ 'Particule Swarm Optimizer (by Leontitsis) [fminswarm]' ];
  options.optimizer = mfilename;
  pars=options;
  return
end
if nargin <= 2
	options=[];
end
if isempty(options)
  options=feval(mfilename, 'defaults');
end
options.Hybrid='none';

[pars,fval,exitflag,output] = fminswarmhybrid(fun, pars, options, varargin{:});

% best optimizers:
% fminimfil
% n<=3 fminsimpsa
% n<=5 fminpso



