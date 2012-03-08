function [pars_out,criteria,message,output] = fits(a, model, pars, options, constraints, varargin)
% [pars,criteria,message,output] = fits(a, model, pars, options, ...) : fit data set on a model
%
%   @iData/fits find best parameters estimates in order to minimize the 
%     fitting criteria using function 'fun' as model, by mean of an optimization method
%     described with the 'options' structure argument.
%     Additional constraints may be set by fxing some parameters, or define
%     more advanced constraints (min, max, steps). The last arguments controls the fitting
%     options with the optimset mechanism, and the constraints to apply during optimization.
%     The fit can be applied sequentially and independently onto iData object arrays.
%  [pars,...] = fits(a, model, pars, options, lb, ub)
%     uses lower and upper bounds as parameter constraints (double arrays)
%  [pars,...] = fits(a, model, pars, options, fixed)
%     indicates which parameters are fixed (non zero elements of array).
%  [pars,...] = fits(a, model, pars, 'optimizer', ...)
%     uses a specific optimizer and its default options.
%  [pars,...] = fits(a, model, pars, options, constraints, args...)
%     send additional arguments to the fit model(pars, axes, args...)
%  [optimizers,functions] = fits(iData)
%     returns the list of all available optimizers and fit functions.
%  fits(iData)
%     displays the list of all available optimizers and fit functions.
%  You may create new fit functions with the 'ifitmakefunc' tool.
% 
% The default fit options.criteria is 'least_square', but others are available:
%   least_square          (|Signal-Model|/Error).^2     non-robust 
%   least_absolute         |Signal-Model|/Error         robust
%   least_median    median(|Signal-Model|/Error)        robust, scalar
%   least_max          max(|Signal-Model|/Error)        non-robust, scalar
%
% input:  a: object or array (iData)
%           when given as an epty iData, the list of optimizers and fit models is shown.
%         model: model function (char/cellstr)
%           when given as a cellstr, the product of all functions is used.
%           when set to empty, the 'gauss' 1D function is used (and possibly extended to multidimensional).
%         pars: initial model parameters (double array). 
%           when set to empty the starting parameters are guessed.
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
%             Optimization method. Default is 'fminsearch' (char/function handle)
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
%           modelValue:        Last model evaluation (iData)
%           parsHistoryUncertainty: Uncertainty on the parameters obtained from 
%                              the optimization trajectory (double)
%
% ex:     p=fits(a,'gauss',[1 2 3 4]);
%         o=fminpowell('defaults'); o.OutputFcn='fminplot'; 
%         [p,c,m,o]=fits(a,'gauss',[1 2 3 4],o); b=o.modelValue; plot(a,b)
%
% Version: $Revision: 1.44 $
% See also iData, fminsearch, optimset, optimget, ifitmakefunc, iFuncs

% private functions: eval_criteria, least_square
 
% handle default parameters, if missing
if nargin == 1
  if isempty(a)
    % return the list of all available optimizers and fit functions
    output = {};
    pars_out   = {};
    if nargout == 0
      fprintf(1, '\n%s\n', version(iData));
      
      fprintf(1, '      OPTIMIZER DESCRIPTION [%s]\n', 'iFit/iOptim');
      fprintf(1, '-----------------------------------------------------------------\n'); 
    end
    d = dir([ fileparts(which(mfilename)) filesep '..' filesep 'iOptim' ]);
    for index=1:length(d)
      this = d(index);
      try
        [dummy, method] = fileparts(this.name);
        options = feval(method,'defaults');
        if isstruct(options)
          output{end+1} = options;
          pars_out{end+1}   = method;
          if nargout == 0
            fprintf(1, '%15s %s\n', options.optimizer, options.algorithm);
          end
        end
      end
    end % for
    if nargout == 0
      fprintf(1, '\n');
      fprintf(1, '       FUNCTION DESCRIPTION [%s]\n', 'iFit/iFuncs');
      fprintf(1, '-----------------------------------------------------------------\n'); 
    end
    d = dir([ fileparts(which(mfilename)) filesep '..' filesep 'iFuncs' ]);
    criteria = {}; 
    for index=1:length(d)
      this = d(index);
      try
        [dummy, method] = fileparts(this.name);
        options = feval(method,'identify');
        if isstruct(options)
          criteria{end+1}   = method;
          if nargout == 0
            fprintf(1, '%15s %s\n', method, options.Name);
          end
        end
      end
    end % for
    % local (pwd) functions
    
    message = '';
    d = dir(pwd);
    for index=1:length(d)
      this = d(index);
      try
        [dummy, method] = fileparts(this.name);
        options = feval(method,'identify');
        if isstruct(options)
          criteria{end+1}   = method;
          if isempty(message)
            fprintf(1, '\nLocal functions in: %s\n', pwd);
            message = ' ';
          end
          if nargout == 0
            fprintf(1, '%15s %s\n', method, options.Name);
          end
        end
      end
    end % for
    if nargout == 0
      fprintf(1, '\n');
      % plot all functions
      m = ceil(sqrt(length(criteria)));
      n = ceil(length(criteria)/m);
      for index=1:length(criteria)
        this = criteria{index};
        subplot(m,n,index); %,'align');
        id   = feval(this, 'plot');
        set(gca,'XTickLabel',[],'XTick',[]); 
        set(gca,'YTickLabel',[],'YTick',[]); 
        xlabel(' '); ylabel(' ');
        axis tight
      end
    end
    message = 'Optimizers and fit functions list'; 
    return
  end
end

if nargin < 2
  model = '';
end
if isempty(model)
  model = 'gauss';
end
if nargin < 3
  pars = [];            % guess starting parameters
end
if nargin < 4, options=[]; end
if isempty(options)
  options = 'fminimfil';% default optimizer
end
if nargin < 5
  constraints = [];     % no constraints
end
if nargin < 6
  varargin = {};
end
if (ischar(options) && length(strtok(options,' =:;'))==length(options)) | isa(options, 'function_handle')
  algo = options;
  options           = feval(algo,'defaults');
  if isa(algo, 'function_handle'), algo = func2str(algo); end
  options.optimizer = algo;
elseif ischar(options), options=str2struct(options);
end

% handle input iData arrays
if numel(a) > 1
  pars_out={}; criteria=[]; message={}; output={};
  for index=1:numel(a)
    [pars_out{index}, criteria(index), message{index}, output{index}] = ...
      fits(a(index), model, pars, constraints, options, varargin{:});
  end
  pars = pars_out;
  return
end

% test the model: is this a function handle or a char ?
if ~isa(model, 'function_handle') && ~ischar(model) && ~iscellstr(model)
  iData_private_error(mfilename,[ 'The model argument is of class ' class(model) '. Should be a function name, expression or function handle.' ]);
end
if ischar(model)
  % is this an expression ?
  if ~exist(model)
    model = ifitmakefunc(model);
    t=which(char(model)); % force to rehash the new function
  end
end

if isempty(model)
  iData_private_error(mfilename,[ 'The model argument is empty. Should be a function name, expression or function handle.' ...
    sprintf('\n') 'Type "fits(iData)" to get a list of available predefined models.' ...
    sprintf('\n') 'or use "ifitmakefunc" to create one.' ]);
end

% handle options
if ~isfield(options, 'optimizer')
  options.optimizer = 'fminpowell'; %  this one does not set 'optimizer'.
end
if ~isfield(options, 'criteria')
  options.criteria  = @least_square;
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
% get starting parameters (guess) if needed
try
  [dummy, info] = ieval(a, model,'guess', varargin{:});    % model info with guessed parameters
catch
  [dummy, info] = ieval(a, model,'identify', varargin{:}); % model info 
end
pars_isstruct=0;
if ischar(pars)
  pars = str2struct(pars);
end
if isstruct(pars)
  % search 'pars' names in the model parameters, and reorder the parameter vector
  for index=1:length(info.Parameters)
    match = strcmp(info.Parameters{index}, fieldnames(pars));
    if ~isempty(match) && any(match)
      p(index) = pars.(info.Parameters{index});
    end
  end
  if length(p) ~= length(info.Parameters)
    disp('Actual parameters')
    disp(pars)
    disp([ 'Required model ' info.Name ' parameters' ])
    disp(info.Parameters)
    iData_private_error(mfilename,[ 'The parameters entered as a structure do not define all required model parameters.' ]);
  else
    pars_isstruct=1;
    pars = p;
  end
elseif isempty(pars)
  pars=info.Guess;               % get default starting parameters
end

if ~isstruct(constraints) && ~isempty(constraints)
  iData_private_error(mfilename,[ 'The constraints argument is of class ' class(constraints) '. Should be a single array or a struct' ]);
end

% removes warnings
iData_private_warning('enter', mfilename);

if ~isfield(options,'Display') options.Display=''; end
if ~isfield(options,'algorithm') options.algorithm=options.optimizer; end

pars = reshape(pars, [ 1 numel(pars)]); % a single row
constraints.parsStart      = pars;
constraints.parsHistory    = [];
constraints.parsNames      = info.Parameters;
constraints.criteriaHistory= [];
constraints.modelName      = info.Name;
constraints.algorithm      = options.algorithm;
constraints.optimizer      = options.optimizer;
constraints.funcCount      = 0;

if isempty(options.Display)    options.Display='notify'; end

if strcmp(options.Display, 'iter') | strcmp(options.Display, 'final') | strcmp(options.Display, 'notify')
  disp([ '** Starting fit of ' a.Tag ' using model ' info.Name ' with optimizer ' options.algorithm ]);
  disp(char(a))
  disp(  '** Minimization performed on parameters:');
  for index=1:length(info.Parameters); fprintf(1,'%10s ', info.Parameters{index}); end; fprintf(1,'\n');
  fprintf(1,'%10.2g ', pars); fprintf(1,'\n');
end

% call minimizer ===============================================================
if abs(nargin(options.optimizer)) == 1 || abs(nargin(options.optimizer)) >= 6
  [pars_out,criteria,message,output] = feval(options.optimizer, ...
    @(pars) eval_criteria(pars, model, options.criteria, a, varargin{:}), pars, options, constraints);
else
  % Constraints not supported by optimizer
  [pars_out,criteria,message,output] = feval(options.optimizer, ...
    @(pars) eval_criteria(pars, model, options.criteria, a, varargin{:}), pars, options);
end

if nargout > 3
  % set output/results
  if ischar(message) | ~isfield(output, 'message')
    output.message = message;
  else
    output.message = [ '(' num2str(message) ') ' output.message ];
  end
  output.modelName  = constraints.modelName;
  output.modelInfo  = info;
  output.modelValue = ieval(a, model, pars_out, varargin{:}); % evaluate model iData
  output.parsNames  = constraints.parsNames;
  output.corrcoef   = eval_corrcoef(a, output.modelValue);
end

if strcmp(options.Display, 'iter') | strcmp(options.Display, 'final') | strcmp(options.Display, 'notify')
  disp([ sprintf('\n') '** Ending fit of ' a.Tag ' using model ' info.Name ' with optimizer ' options.algorithm ]);
  fprintf(1, '   %i iterations. Criteria=%g Status: %s\n', output.funcCount, sum(criteria(:)), output.message);
  for index=1:length(info.Parameters); fprintf(1,'%10s ', info.Parameters{index}); end; fprintf(1,'\n');
  fprintf(1,'%10.2g ', pars_out); fprintf(1,'\n');
  index=find(output.criteriaHistory < min(output.criteriaHistory)*4);   % identify tolerance region around optimum
  
  % write uncertainty if not done by fmin_private_wrapper
  if ~strcmp(options.Display,'final') & ~strcmp(options.Display,'iter')
    if length(index) > 10 % test length of tolerence region
      disp('** Gaussian uncertainty on parameters (half width, from the optimization history)')
      fprintf(1,'%10.2g ', output.parsHistoryUncertainty); fprintf(1,'\n');
    end
    if isfield(output,'parsHessianUncertainty') && ~isempty(output.parsHessianUncertainty)
      fprintf(1,'** Gaussian uncertainty on parameters (half width, from the Hessian matrix)');
      if isfield(output, 'corrcoef')
        fprintf(1, ', CorrCoef=%g\n', output.corrcoef);
      else disp(' '); end
      fprintf(1,'%10.2g ', output.parsHessianUncertainty); fprintf(1,'\n');
    end
  end
end

% reset warnings
iData_private_warning('exit', mfilename);

end % fits end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = eval_criteria(pars, model, criteria, a, varargin)
% criteria to minimize
  if nargin<5, varargin={}; end
  % then get data
  Signal = iData_private_cleannaninf(get(a,'Signal'));
  Error  = iData_private_cleannaninf(get(a,'Error'));
  Model  = ieval(a, model, pars, varargin{:}); % return signal=model values*monitor and monitor
  Model  = iData_private_cleannaninf(get(Model, 'Signal'));
  if isempty(Model)
    iData_private_error(mfilename,[ 'The model ' model ' could not be evaluated (returned empty).' ]);
  end
  m      = iData_private_cleannaninf(get(a,'Monitor')); m=real(m);
  if not(all(m == 1 | m == 0)),
    Model  = genop(@rdivide,Model,m);            % fit(signal/monitor) 
    Signal = genop(@rdivide,Signal,m); Error=genop(@rdivide,Error,m); % per monitor
  end
  
  % compute criteria
  c = feval(criteria, Signal(:), Error(:), Model(:));
  % devide by the number of degrees of freedom
  % <http://en.wikipedia.org/wiki/Goodness_of_fit>
  c = c/(numel(Signal) - length(pars) - 1); % reduced 'Chi^2'
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
function r=eval_corrcoef(a, modelValue)
% correlation coefficient between the data and the model

  Signal = iData_private_cleannaninf(get(a,'Signal'));
  Error  = iData_private_cleannaninf(get(a,'Error'));
  Model  = iData_private_cleannaninf(get(modelValue, 'Signal'));
  m      = iData_private_cleannaninf(get(a,'Monitor')); m=real(m);
  if not(all(m(:) == 1 | m(:) == 0)),
    Model  = genop(@rdivide,Model,m);            % fit(signal/monitor) 
    Signal = genop(@rdivide,Signal,m); Error=genop(@rdivide,Error,m); % per monitor
  end
  
  % compute the correlation coefficient
  if isempty(Error) || isscalar(Error) || all(Error(:) == Error(end))
    wt = 1;
  else
    wt = 1./Error;
    wt(find(~isfinite(wt))) = 0;
  end
  r  = corrcoef(Signal.*wt,Model.*wt);
  r  = r(1,2);                                     % correlation coefficient
end % eval_corrcoef

