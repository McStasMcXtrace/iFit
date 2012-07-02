function [pars_out,criteria,message,output] = fits(model, a, pars, options, constraints, varargin)
% [pars,criteria,message,output] = fits(model, data, pars, options, ...) : fit a model on a data set
%
%   @iFunc/fits find best parameter estimates in order to minimize the 
%     fitting criteria using model 'fun', by mean of an optimization method
%     described with the 'options' structure argument.
%     Additional constraints may be set by fixing some parameters, or define
%     more advanced constraints (min, max, steps). The last arguments controls the fitting
%     options with the optimset mechanism, and the constraints to apply during optimization.
%  [pars,...] = fits(a, model, pars, options, lb, ub)
%     uses lower and upper bounds as parameter constraints (double arrays)
%  [pars,...] = fits(a, model, pars, options, fixed)
%     indicates which parameters are fixed (non zero elements of array).
%  [pars,...] = fits(a, model, pars, 'optimizer', ...)
%     uses a specific optimizer and its default options options=feval(optimizer,'defaults')
%  [pars,...] = fits(a, model, pars, options, constraints, args...)
%     send additional arguments to the fit model(pars, axes, args...)
% 
% The default fit options.criteria is 'least_square', but others are available:
%   least_square          (|Signal-Model|/Error).^2     non-robust 
%   least_absolute         |Signal-Model|/Error         robust
%   least_median    median(|Signal-Model|/Error)        robust, scalar
%   least_max          max(|Signal-Model|/Error)        non-robust, scalar
%
%  Type <a href="matlab:doc(iData,'Fit')">doc(iData,'Fit')</a> to access the iFit/Fit Documentation.
%
% input:  model: model function (iFunc)
%         data: array or structure/object (numeric or structure or cell)
%               Can be entered as a single numeric array (the Signal), or as a structure/object
%                 with possible members Signal, Error, Monitor, Axes={x,y,...}
%               or as a cell { x,y, ... , Signal }
%         pars: initial model parameters (double array). 
%           when set to empty or 'guess', the starting parameters are guessed.
%         options: structure as defined by optimset/optimget (char/struct)
%           if given as a char, it defines the algorithm to use and its default options.
%           when set to empty, it sets the default algorithm options (fminimfil).
%           options.TolX
%             The termination tolerance for x. Its default value is 1.e-4.
%           options.TolFun
%             The termination tolerance for the function value. The default value is 1.e-4. 
%             This parameter is used by fminsearch, but not fminbnd.
%           options.MaxIter
%             Maximum number of iterations allowed.
%           options.MaxFunEvals
%             The maximum number of function evaluations allowed. 
%           options.optimizer
%             Optimization method. Default is 'fminpowell' (char/function handle)
%             the syntax for calling the optimizer is e.g. optimizer(criteria,pars,options,constraints)
%           options.criteria
%             Minimization criteria. Default is 'least_square' (char/function handle)
%             the syntax for evaluating the criteria is criteria(Signal, Error, Model)
%           options.OutputFcn
%             Function called at each iteration as outfun(pars, optimValues, state)
%             The 'fminplot' function may be used.
%           options.Display
%             Display additional information during fit: 'iter','off','final'. Default is 'iter'.
%         constraints: fixed parameter array. Use 1 for fixed parameters, 0 otherwise (double array)
%           OR use empty to not set constraints
%           OR use a structure with some of the following fields:
%           constraints.min:   minimum parameter values (double array)
%           constraints.max:   maximum parameter values (double array)
%           constraints.step:  maximum parameter step/change allowed.
%           constraints.fixed: fixed parameter flag. Use 1 for fixed parameters, 0 otherwise (double array)
%
% output: 
%         pars:              best parameter estimates (double array)
%         criteria:          minimal criteria value achieved (double)
%         message:           return message/exitcode from the optimizer (char/integer)
%         output:            additional information about the optimization (structure)
%           algorithm:         Algorithm used (char)
%           funcCount:         Number of function evaluations (double)
%           iterations:        Number of iterations (double)
%           parsHistory:       Parameter set history during optimization (double array)
%           criteriaHistory:   Criteria history during optimization (double array)
%           modelValue:        Final best model evaluation
%           parsHistoryUncertainty: Uncertainty on the parameters obtained from 
%                              the optimization trajectory (double)
%
% ex:     p=fits(gauss, data,[1 2 3 4]);
%         o=fminpowell('defaults'); o.OutputFcn='fminplot'; 
%         [p,c,m,o]=fits(gauss,data,[1 2 3 4],o); b=o.modelValue; plot(a,b)
%
% Version: $Revision: 1.1 $
% See also fminsearch, optimset, optimget, iFunc

% first get the axes and signal from 'data'

% a.Signal (numeric)
% a.Error (numeric)
% a.Monitor (numeric)
% a.Axes (cell of numeric)

% check of input arguments =====================================================

if isempty(model)
  disp([ 'iFunc:' mfilename ': Using default gaussian model as fit function.' ]);
  model = gauss;
end

if nargin < 2
	a = [];
end

% extract Signal from input argument, as well as a Data identifier
% default values
Monitor=1; Error=1; Axes={}; Signal=[]; Name = '';
if iscell(a)
  Signal = a{end};
  a(end) = [];
  Axes = a;
elseif isstruct(a) || isa(a, 'iData')
  if isfield(a,'Signal')  Signal  = a.Signal; end
  if isfield(a,'Error')   Error   = a.Error; end
  if isfield(a,'Monitor') Monitor = a.Monitor; end
  if isfield(a,'Axes')    Axes    = a.Axes; end
  if isa(a, 'iData')
    Axes=cell(1,ndims(a));
    for index=1:ndims(a)
      Axes{index} = getaxis(a, index);
    end
    Name = [ inputname(2) ' ' char(a) ];
  end
elseif isnumeric(a)
  Signal = a; 
end
if isempty(Name)
  Name   = [ class(a) ' ' mat2str(size(Signal)) ' ' inputname(2) ];
end

% handle parameters: from char, structure or vector
if nargin < 3
  pars = []; % will use default/guessed parameters
end
pars_isstruct=0;
if ischar(pars) && ~strcmp(pars,'guess')
  pars = str2struct(pars);
end
if isstruct(pars)
  % search 'pars' names in the model parameters, and reorder the parameter vector
  p = [];
  for index=1:length(model.Parameters)
    match = strcmp(model.Parameters{index}, fieldnames(pars));
    if ~isempty(match) && any(match)
      p(index) = pars.(model.Parameters{index});
    end
  end
  if length(p) ~= length(model.Parameters)
    disp('Actual parameters')
    disp(pars)
    disp([ 'Required model ' model.Name ' ' model.Tag ' parameters' ])
    disp(model.Parameters)
    error([ 'iFunc:' mfilename], [ 'The parameters entered as a structure do not define all required model parameters.' ]);
  else
    pars_isstruct=1;
    pars = p;
  end
elseif strcmp(pars,'guess')
  feval(model, pars, Axes{:}, Signal);           % get default starting parameters
  pars = model.ParameterValues;
elseif isempty(pars)
  if ~isempty(model.ParameterValues)
    pars = model.ParameterValues;     % use stored starting parameters
  else     
    feval(model, pars, Axes{:}, Signal);         % get default starting parameters
    pars = model.ParameterValues;
  end
end
pars = reshape(pars, [ 1 numel(pars)]); % a single row

% handle options
if nargin < 4, options=[]; end
if isempty(options)
  options = 'fminpowell';% default optimizer
end
if (ischar(options) && length(strtok(options,' =:;'))==length(options)) | isa(options, 'function_handle')
  algo = options;
  options           = feval(algo,'defaults');
  if isa(algo, 'function_handle'), algo = func2str(algo); end
  options.optimizer = algo;
elseif ischar(options), options=str2struct(options);
end
if ~isfield(options, 'optimizer')
  options.optimizer = 'fminpowell';
end
if ~isfield(options, 'criteria')
  options.criteria  = @least_square;
end
if ~isfield(options,'Display')   options.Display=''; end
if isempty(options.Display)    options.Display='notify'; end
if ~isfield(options,'algorithm') options.algorithm=options.optimizer; end

% handle constraints
if nargin < 5
  constraints = [];     % no constraints
end
% handle constraints given as vectors
if (length(constraints)==length(pars) | isempty(pars)) & (isnumeric(constraints) | islogical(constraints))
  if nargin<6
    fixed            = constraints;
    constraints      =[];
    constraints.fixed=fixed;
  elseif isnumeric(varargin{1}) && ~isempty(varargin{1}) ...
      && length(constraints) == length(varargin{1})
    % given as lb,ub parameters (nargin==6)
    lb = constraints;
    ub = varargin{1};
    varargin(1) = []; % remove the 'ub' from the additional arguments list
    constraints     = [];
    constraints.min = lb;
    constraints.max = ub;
  end
end
if ~isstruct(constraints) && ~isempty(constraints)
  error([ 'iFunc:' mfilename],[ 'The constraints argument is of class ' class(constraints) '. Should be a single array or a struct' ]);
end
constraints.parsStart      = pars;
constraints.parsHistory    = [];
constraints.criteriaHistory= [];
constraints.algorithm      = options.algorithm;
constraints.optimizer      = options.optimizer;
constraints.funcCount      = 0;

% handle arrays of model functions
if numel(model) > 1
  pars_out={} ; criteria={}; message={}; output={};
  for index=1:numel(model)
    [pars_out{end+1},criteria{end+1},message{end+1},output{end+1}]=fits(model(index), a, pars, options, constraints, varargin{:});
  end
  return
end

if isempty(Signal)
  error([ 'iFunc:' mfilename ],[ 'Undefined/empty Signal ' inputname(2) ' to fit. Syntax is fits(model, Signal, parameters, ...).' ]);
end

if isvector(Signal) ndimS = 1;
else                ndimS = ndims(Signal)
end
% handle case when model dimensionality is larger than actual Signal
if model.Dimension > ndimS
  error([ 'iFunc:' mfilename ], 'Signal %s with dimensionality %d has lower dimension than model %s dimensionality %d\n', Name, ndimS, model.Name, model.Dimension);
% handle case when model dimensionality is smaller than actual Signal
elseif model.Dimension < ndimS && rem(ndimS, model.Dimension) == 0
  % extend model to match Signal dimensions
  disp(sprintf('iFunc:%s: Extending model %s dimensionality %d to data %s dimensionality %d.\n', ...
    mfilename, model.Name, model.Dimension, Name, ndimS));
  new_model=model;
  for index=2:(ndimS/model.Dimension)
    new_model = new_model * model;
  end
  model = new_model;
  clear new_model
elseif model.Dimension ~= ndimS
  error([ 'iFunc:' mfilename ], 'Signal %s with dimensionality %d has higher dimension than model %s dimensionality %d.\n', Name, ndimS, model.Name, model.Dimension);
end

% create the new Data structure to pass to the criteria
a = [];
a.Signal = iFunc_private_cleannaninf(Signal);
a.Error  = iFunc_private_cleannaninf(Error);
a.Monitor= iFunc_private_cleannaninf(Monitor);
a.Axes   = Axes;
clear Signal Error Monitor Axes

% starting configuration
feval(model, pars, a.Axes{:}, a.Signal); % this updates the 'model' with starting parameter values

% we need to call the optimization method with the eval_criteria as FUN
% call minimizer ===============================================================
if abs(nargin(options.optimizer)) == 1 || abs(nargin(options.optimizer)) >= 6
  [pars_out,criteria,message,output] = feval(options.optimizer, ...
    @(pars) eval_criteria(model, pars, options.criteria, a, varargin{:}), pars, options, constraints);
else
  % Constraints not supported by optimizer
  [pars_out,criteria,message,output] = feval(options.optimizer, ...
    @(pars) eval_criteria(model, pars, options.criteria, a, varargin{:}), pars, options);
end

% format output arguments ======================================================
pars_out = reshape(pars_out, [ 1 numel(pars_out) ]); % row vector
model.ParameterValues = pars_out;
if nargout > 3
  output.modelValue = feval(model, pars_out, a.Axes{:});
  output.modelAxes  = a.Axes{:};
  output.model      = model;
  
  % set output/results
  if ischar(message) | ~isfield(output, 'message')
    output.message = message;
  else
    output.message = [ '(' num2str(message) ') ' output.message ];
  end
  output.parsNames  = model.Parameters;
  output.corrcoef   = eval_corrcoef(a.Signal, a.Error, a.Monitor, output.modelValue);
end
if pars_isstruct
  pars_out = cell2struct(num2cell(pars_out), strtok(model.Parameters), 2);
end

end % fits




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = eval_criteria(model, pars, criteria, a, varargin)
% criteria to minimize
  if nargin<5, varargin={}; end
  % then get model value
  Model  = feval(model, pars, a.Axes{:}, varargin{:}); % return model values
  Model  = iFunc_private_cleannaninf(Model);
  if isempty(Model)
    error([ 'iFunc:' mfilename ],[ 'The model ' model ' could not be evaluated (returned empty).' ]);
  end
  a.Monitor =real(a.Monitor);
  if not(all(a.Monitor == 1 | a.Monitor == 0)),
    Model    = bsxfun(@rdivide,Model,   a.Monitor); % fit(signal/monitor) 
    a.Signal = bsxfun(@rdivide,a.Signal,a.Monitor); 
    a.Error  = bsxfun(@rdivide,a.Error, a.Monitor); % per monitor
  end
  
  % compute criteria
  c = feval(criteria, a.Signal(:), a.Error(:), Model(:));
  % divide by the number of degrees of freedom
  % <http://en.wikipedia.org/wiki/Goodness_of_fit>
  if numel(a.Signal) > length(pars)-1
    c = c/(numel(a.Signal) - length(pars) - 1); % reduced 'Chi^2'
  end
end % eval_criteria

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <http://en.wikipedia.org/wiki/Least_squares>
function c=least_square(Signal, Error, Model)
% weighted least square criteria, which is also the Chi square
% the return value is a vector, and most optimizers use its sum (except LM).
% (|Signal-Model|/Error).^2
  c = least_absolute(Signal, Error, Model);
  c = c.*c;
end % least_square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <http://en.wikipedia.org/wiki/Least_absolute_deviation>
function c=least_absolute(Signal, Error, Model)
% weighted least absolute criteria
% the return value is a vector, and most optimizers use its sum (except LM).
% |Signal-Model|/Error
  if isempty(Error) || isscalar(Error) || all(Error == Error(end))
    index = find(isfinite(Model) & isfinite(Signal));
    c = abs(Signal(index)-Model(index)); % raw least absolute
  else
    % find minimal non zero Error
    Error = abs(Error);
    index = find(Error~=0 & isfinite(Error));
    minError = min(Error(index));
    % find zero Error, which should be replaced by minimal Error whenever possible
    index = find(Error == 0);
    Error(index) = minError;
    index = find(isfinite(Error) & isfinite(Model) & isfinite(Signal));
    if isempty(index), c=Inf;
    else               c=abs((Signal(index)-Model(index))./Error(index));
    end
  end
end % least_absolute

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <http://en.wikipedia.org/wiki/Median_absolute_deviation>
function c=least_median(Signal, Error, Model)
% weighted median absolute criteria
% the return value is a scalar
% median(|Signal-Model|/Error)
  c = median(least_absolute(Signal, Error, Model));
end % least_median

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <http://en.wikipedia.org/wiki/Absolute_deviation>
function c=least_max(Signal, Error, Model)
% weighted median absolute criteria
% the return value is a scalar
% median(|Signal-Model|/Error)
  c = max(least_absolute(Signal, Error, Model));
end % least_max

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=eval_corrcoef(Signal, Error, Monitor, Model)
% correlation coefficient between the data and the model

  if not(all(Monitor(:) == 1 | Monitor(:) == 0)),
    Model  = bsxfun(@rdivide,Model, Monitor); % fit(signal/monitor) 
    Signal = bsxfun(@rdivide,Signal,Monitor); 
    Error  = bsxfun(@rdivide,Error, Monitor); % per monitor
  end
  
  % compute the correlation coefficient
  if isempty(Error) || isscalar(Error) || all(Error(:) == Error(end))
    wt = 1;
  else
    wt = 1./Error;
    wt(find(~isfinite(wt))) = 0;
  end
  r = corrcoef(Signal.*wt,Model.*wt);
  r = r(1,2);                                     % correlation coefficient
  if isnan(r)
    r = corrcoef(Signal,Model);
    r = r(1,2);
  end
end % eval_corrcoef

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=iFunc_private_cleannaninf(s)
% iFunc_private_cleannaninf: clean NaNs and Infs from a numerical field
%
  
  if isnumeric(s)
    S = s(:);
    if all(isfinite(S)), return; end
    index_ok     = find(isfinite(S));

    maxs = max(S(index_ok));
    mins = min(S(index_ok));

    S(isnan(S)) = 0;
    if ~isempty(mins)
      if mins<0, S(find(S == -Inf)) = mins*100;
      else       S(find(S == -Inf)) = mins/100; end
    end
    if ~isempty(maxs)
      if maxs>0, S(find(S == +Inf)) = maxs*100;
      else       S(find(S == +Inf)) = maxs/100; end
    end

    s = double(reshape(S, size(s)));
  end

end
