function [pars,fval,exitflag,output] = fmax(objective, pars, options,  varargin)
% [Maxmimum,FVAL,EXITFLAG,OUTPUT] = fmax(FUN,PARS,[OPTIONS],[X,Y,...]) Best maximizer
%
% This maximization method is determined automatically from the objective function
% behaviour and number of free parameters. You can however force a specific 
% optimizer by setting e.g. options.optimizer='fminpso'
%
% WARNING: as the selected optimizer may change from one call to an other, the
% solution found may vary as well. To avoid that, rather use a specific optimizer.
%
% Syntax:
%   fmax(model)
%     maximise the model, using current model or guessed parameters. Return best parameters
%   fmax(model, pars)
%     maximise the model, starting with given parameters and return best parameter set 
%     when pars is given as [], current/guessed parameters are used
%   fmax(model, pars, 'optimizer')
%     maximise the model, as above, using given optimizer with its default configuration
%     for instance: fmax(model, [], 'fminpso')
%   fmax(model, pars, options)
%     maximise the model, as above, using given optimizer configuration
%     for instance: fmax(model, [], 'optimizer=fminpso; OutputFcn=fminplot; Display=iter')
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
%   fmax(fun, pars) asks to maximize the 'fun' iFunc model with starting
%     parameters 'pars' (vector)
%   fmax(fun, pars, 'optimizer') 
%     use optimizer with its default options
%   fmax(fun, pars, options) 
%     same as above, with customized options (optimset)
%   fmax(fun, pars, options, x,y,...) 
%     same as above, with specific axes
%
% The options structure may contain the following members, in agreement with 'optimset':
%    options.Display: Level of display [ off | iter | notify | final ]. Default is 'off'
%    options.MaxFunEvals: Maximum number of function evaluations allowed, sometimes 
%      referred as the 'cost' or 'budget'.
%    options.MaxIter: Maximum number of iterations allowed
%    options.TolFun: Termination tolerance on the function value (absolute value or change). 
%      Use 'x%' to specify a relative function change.
%    options.TolX: Termination tolerance on parameter change. 
%      Use 'x%' to specify a relative parameter change.
%    options.OutputFcn: Name of an output function. When set, it is called at each
%      iteration step. You may use 'fminplot', which is provided in Optimizers. 
%      Refer to the Fit page for more information about fminplot. A simpler/faster
%      alternative is the 'fminstop' option.
%    options.PlotFcns: same as OutputFcn, but can be a set of function in a cell array.
%    options.FunValCheck: Check for invalid values, such as NaN or complex
%    options.MinFunEvals: when set, waits for a given number of iterations before testing for convergence
%    options.optimizer: the optimizer to use
%
% Input:
%  FUN is the iFunc model to maximize (handle or string)
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess. PARS can also be given as a single-level structure.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%     o=fmax('defaults')
%  options.MinFunEvals sets the minimum number of function evaluations to reach
%  An empty OPTIONS sets the default configuration. 
%
%  X,Y,...: additional axes to use.
%
% Output:
%          MAXIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MAXIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% Example:
%  model= gauss1;
%  fix(model, 'all'); model.Intensity='free';
%  model.Intensity=1; model.HalfWidth=.5;
%  xlim(model, 'Intensity',[-2 2])
%  fmax(model)
%
% Version: $Date$
% See also: fminsearch, optimset, iFunc/fmin
% (c) E.Farhi, ILL. License: EUPL.

% we maximize the iFunc: (p)feval(iFunc, p). Must guess some axes to use.


objective2 = -objective;

if nargin < 2, pars    = []; end
if nargin < 3, options = []; end

[pars,fval,exitflag,output] = fmin(objective2, pars, options,  varargin{:});

if ~isempty(inputname(1))
  objective.UserData        = objective2.UserData;
  objective.ParameterValues = objective2.ParameterValues;
  objective.Constraint      = output.constraints; % restore initial constraints
  assignin('caller', inputname(1), objective);
end

