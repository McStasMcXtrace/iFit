function [pars,fval,exitflag,output] = fmin(objective, pars, options,  varargin)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = FMIN(FUN,PARS,[OPTIONS],[X,Y,...]) Best optimizer
%
% This minimization method is determined automatically from the objective function
% behaviour and number of free parameters. You can however force a specific 
% optimizer by setting e.g. options.optimizer='fminpso'
%
% Syntax:
%   fmin(fun, pars) asks to minimize the 'fun' iFunc model with starting
%     parameters 'pars' (vector)
%   fmin(fun, pars, 'optimizer') 
%     use optimizer with its default options, 
%     e.g: fmin(fun, [], 'fminpso')
%   fmin(fun, pars, options) 
%     same as above, with customized options (optimset), 
%     e.g.: fmin(fun, [], 'optimizer=fminpso; OutputFcn=fminplot; Display=iter')
%   fmin(fun, pars, options, x,y,...) 
%     same as above, with specific axes
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
%  FUN is the iFunc model to minimize (handle or string)
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
%  X,Y,...: additional axes to use.
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% Example:
%  model= gauss1;
%  fix(model, 'all'); model.Intensity='free';
%  model.Intensity=1; model.HalfWidth=.5;
%  xlim(model, 'Intensity',[-2 2])
%  fmin(model)
%
% Version: $Date$
% See also: fminsearch, optimset, iFunc/fmax
% (c) E.Farhi, ILL. License: EUPL.

% we minimize the iFunc: (p)feval(iFunc, p). Must guess some axes to use.

if nargin < 2, pars    = []; end
if nargin < 3, options = []; end

% handle parameters: from char, structure or vector
[pars, pars_isstruct, ax] = iFunc_private_get_pars(objective, pars, varargin);

% create the evaluation of the function
fun = @(p)feval(objective', p, ax{:});

% carry 'constraints' from the iFunc if not in input

constraints = objective.Constraint;
objective.Constraint=[];  % we have transfered the restraints.

[pars,fval,exitflag,output] = fmin(fun, pars, options, constraints);

if ~isempty(inputname(1))
  objective.UserData.output = output;
  objective.ParameterValues = pars(:);
  objective.Constraint      = constraints; % restore initial constraints
  assignin('caller', inputname(1), objective);
end
% return struct when pars where given as such,only those changed
if ~isempty(pars_isstruct)
  pars = pars(pars_isstruct);
  pars_name = objective.Parameters(pars_isstruct);
  pars = cell2struct(num2cell(pars(:)), strtok(pars_name(:)), 1);
end
