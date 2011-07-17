function [pars_out,criteria,message,output] = fits(a, model, pars, options, constraints,ub)
% [pars,criteria,message,output] = fits(a, model, pars, options, constraints) : fit data set on a model
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
%  [optimizers,functions] = fits(iData)
%     returns the list of all available optimizers and fit functions.
%  fits(iData)
%     displays the list of all available optimizers and fit functions.
%  You may create new fit functions with the 'ifitmakefunc' tool.
%
% input:  a: object or array (iData)
%         model: model function (char/cellstr)
%           when given as a cellstr, the product of all functions is used
%         pars: initial model parameters (double array)
%         options: structure as defined by optimset/optimget (char/struct)
%           if given as a char, it defines the algorithm to use and its default options
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
%         o=fminimfil('defaults'); o.OutputFcn='fminplot'; 
%         [p,c,m,o]=fits(a,'gauss',[1 2 3 4],o); b=o.modelValue
%
% Version: $Revision: 1.28 $
% See also iData, fminsearch, optimset, optimget, ifitmakefunc

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
        subplot(m,n,index);
        id   = feval(this, 'plot');
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
  pars = [];
end
if nargin < 5
  constraints = [];
end
if nargin < 4, options=[]; end
if isempty(options)
  options = 'fminimfil';
end
if ischar(options) | isa(options, 'function_handle')
  algo = options;
  options           = feval(algo,'defaults');
  options.optimizer = algo;
end

% handle input iData arrays
if length(a) > 1
  pars_out={}; criteria=[]; message={}; output={};
  for index=1:length(a(:))
    [pars_out{index}, criteria(index), message{index}, output{index}] = ...
      fits(a(index), model, pars, constraints, options);
  end
  pars = pars_out;
  return
end

% test the model: is this a function handle or a char ?
if ~isa(model, 'function_handle') && ~ischar(model)
  iData_private_error(mfilename,[ 'The model argument is of class ' class(model) '. Should be a function name, expression or function handle.' ]);
end
if ischar(model)
  % is this an expression ?
  if ~exist(model)
    model = ifitmakefunc(model);
    t=which(char(model));
  end
end

if isempty(model)
  iData_private_error(mfilename,[ 'The model argument is empty. Should be a function name, expression or function handle.' ]);
end

% handle options
if ~isfield(options, 'optimizer')
  options.optimizer = 'fminsearch';
end
if ~isfield(options, 'criteria')
  options.criteria  = @least_square;
end
try
  [dummy, info] = ieval(a, model,'guess');    % model info with guessed parameters
catch
  [dummy, info] = ieval(a, model,'identify'); % model info 
end
if isempty(pars)
  pars=info.Guess;               % get default starting parameters
end
if isnumeric(constraints) | islogical(constraints)
  if nargin<6
    fixed            = constraints;
    constraints      =[];
    constraints.fixed=fixed;
  else
    lb = constraints;
    constraints     = [];
    constraints.min = lb;
    constraints.max = ub;
  end
end
if ~isstruct(constraints)
  iData_private_error(mfilename,[ 'The constraints argument is of class ' class(constraints) '. Should be a single array or a struct' ]);
end

% removes warnings
try
  warn.set = warning('off','iData:setaxis');
  warn.get = warning('off','iData:getaxis');
catch
  warn = warning('off');
end

if ~isfield(options,'Display') options.Display=''; end
if ~isfield(options,'algorithm') options.algorithm=options.optimizer; end
if isempty(options.Display)    options.Display=''; end

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
try
  [pars_out,criteria,message,output] = feval(options.optimizer, ...
    @(pars) eval_criteria(pars, model, options.criteria, a), pars, options, constraints);
catch
  if strcmp(options.Display, 'iter') | strcmp(options.Display, 'final') | strcmp(options.Display, 'notify')
    disp([ '** ' mfilename ': Error or Constraints not supported by optimizer ' options.optimizer ])
  end
  [pars_out,criteria,message,output] = feval(options.optimizer, ...
    @(pars) eval_criteria(pars, model, options.criteria, a), pars, options);
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
  output.modelValue = ieval(a, model, pars_out); % evaluate model iData
  output.parsNames  = constraints.parsNames;
  output.corrcoeff  = eval_corrcoef(a, output.modelValue);
end

if strcmp(options.Display, 'iter') | strcmp(options.Display, 'final') | strcmp(options.Display, 'notify')
  disp([ sprintf('\n') '** Ending fit of ' a.Tag ' using model ' info.Name ' with optimizer ' options.algorithm ]);
  fprintf(1, '   %i iterations. Status: %s\n', output.funcCount, output.message);
  for index=1:length(info.Parameters); fprintf(1,'%10s ', info.Parameters{index}); end; fprintf(1,'\n');
  fprintf(1,'%10.2g ', pars_out); fprintf(1,'\n');
  index=find(output.criteriaHistory < min(output.criteriaHistory)*4);   % identify tolerance region around optimum
  if length(index) > 10
    disp('** Gaussian uncertainty on parameters (half width, from the optimization history)')
    fprintf(1,'%10.2g ', output.parsHistoryUncertainty); fprintf(1,'\n');
  else
    disp('** Gaussian uncertainty on parameters (half width, from the Hessian matrix)')
    fprintf(1,'%10.2g ', output.parsHessianUncertainty); fprintf(1,'\n');
  end
end


% reset warnings
try
  warning(warn.set);
  warning(warn.get);
catch
  warning(warn);
end

end % fits end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = eval_criteria(pars, model, criteria, a)
% criteria to minimize
  
  % then get data
  Signal = iData_private_cleannaninf(get(a,'Signal'));
  Error  = iData_private_cleannaninf(get(a,'Error'));
  Model  = ieval(a, model, pars); % return signal=model values*monitor and monitor
  Model  = iData_private_cleannaninf(get(Model, 'Signal'));
  m      = iData_private_cleannaninf(get(a,'Monitor')); m=real(m);
  if not(all(m == 1 | m == 0)),
    Model  = genop(@rdivide,Model,m);            % fit(signal/monitor) 
    Signal = genop(@rdivide,Signal,m); Error=genop(@rdivide,Error,m); % per monitor
  end
  
  % compute criteria
  c = feval(criteria, Signal(:), Error(:), Model(:));
  c = c/(prod(size(Signal)) - length(pars) - 1); % normalized 'Chi^2'
end % eval_criteria

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c=least_square(Signal, Error, Model)
% weighted least square criteria, which is also the Chi square
  if all(Error == 0)
    c = sum(abs(Signal-Model).^2); % raw least square
  else
    index = find(Error~=0 & ~isnan(Error) & ~isinf(Error));
    c=(Signal(index)-Model(index))./Error(index);
    c=abs(c);
    c=abs(c.*c);                % Chi square
  end
end % least_square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=eval_corrcoef(a, modelValue)
% correlation coefficient between the data and the model

  Signal = iData_private_cleannaninf(get(a,'Signal'));
  Error  = iData_private_cleannaninf(get(a,'Error'));
  Model  = iData_private_cleannaninf(get(modelValue, 'Signal'));
  m      = iData_private_cleannaninf(get(a,'Monitor')); m=real(m);
  if not(all(m == 1 | m == 0)),
    Model  = genop(@rdivide,Model,m);            % fit(signal/monitor) 
    Signal = genop(@rdivide,Signal,m); Error=genop(@rdivide,Error,m); % per monitor
  end
  
  % compute the correlation coefficient
  if all(Error == 0)
    wt = 1;
  else
    wt = 1./Error;
    wt(find(~isfinite(wt))) = 0;
  end
  r  = corrcoef(Signal.*wt,Model.*wt);
  r  = r(1,2);                                     % correlation coefficient
end % eval_corrcoef

