function [pars,fval,exitflag,output] = fminswarm(fun, pars, options, varargin)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = FMINSWARM(FUN,PARS,[OPTIONS],[CONSTRAINTS]) Particle Swarm Optimization
%
% This minimization method uses a Particle Swarm Optimization algorithm for 
% finding the minimum of the function 'FUN' in the real space. 
%
% Calling:
%   fminswarm(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminswarm(fun, pars, options) same as above, with customized options (optimset)
%   fminswarm(fun, pars, options, fixed) 
%     is used to fix some of the parameters. The 'fixed' vector is then 0 for
%     free parameters, and 1 otherwise.
%   fminswarm(fun, pars, options, lb, ub) 
%     is used to set the minimal and maximal parameter bounds, as vectors.
%   fminswarm(fun, pars, options, constraints) 
%     where constraints is a structure (see below).
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminswarm(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%      o=fminswarm('defaults');
%   option.PopulationSize sets the number of particules in the swarm (20-40).
%   option.SwarmC1 sets the local attractors strength (1-3)
%   option.SwarmC2 sets the global attractor strength (1-3).
%   option.SwarmW  sets inertia weight (0-1).
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
% Reference:
% Kennedy J., Eberhart R.C. (1995): Particle swarm optimization. In: Proc.
% IEEE Conf. on Neural Networks, IV, Piscataway, NJ, pp. 1942-1948
% Shi, Y. and Eberhart, R. C. A modified particle swarm optimizer. Proc. 
% IEEE Int Conf. on Evol Comput pp. 69-73. IEEE Press, Piscataway, NJ, 1998
%
% Contrib:
% Alexandros Leontitsis leoaleq@yahoo.com Ioannina, Greece 2004
%
% Version: $Revision: 1.7 $
% See also: fminsearch, optimset

% this is a wrapper to fminswarmhybrid, without hybrid optimizer

if nargin == 1 & strcmp(fun,'defaults')
  options=fminswarmhybrid('defaults');
  options.Hybrid='none';
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


