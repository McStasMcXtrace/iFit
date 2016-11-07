function [signal, model, ax, name] = feval(model, p, varargin)
% [signal, model, axes, name] = feval(model, parameters, x,y, ...) evaluate a function
%
%   @iFunc/feval applies the function 'model' using the specified parameters and axes
%     and function parameters 'pars' with optional additional parameters.
%     a fast notation is to pass arguments directly to the model:
%       model(p, x,y,z,...)
%
%   signal = feval(model, 'guess', x,y, ..., signal...)
%     makes a quick parameter guess. This usually requires to specify the signal
%     to guess from to be passed after the axes. The parameters are then in model.ParameterValues
%     or returned instead of the signal when the model can not be updated.
%   signal = feval(model, NaN, x,y, ..., signal...)
%     same as above, but force to get the evaluated function value with
%     guessed parameters.
%   signal = feval(model, [ ... Nan ... ], x,y, ..., signal...)
%     requires some of the initial parameters to be given, others as NaN's. These
%     values are then replaced by guessed ones, and the model value is returned.
%   signal = feval(model, parameters, x,y, ...)
%     evaluates the model with given parameters and axes
%   signal = model(p, iData object)
%     evaluates the model on the given iData object axes
%
% input:  model: model function (iFunc, single or array)
%         parameters: model parameters (vector, cell or vectors, structure, iData) or 'guess'
%         x,y,..:  axes values to be used for the computation (vector,matrix,iData)
%         ...: additional parameters may be passed, which are then forwarded to the model
% output: signal: result of the evaluation (vector/matrix/cell) or guessed parameters (vector)
%         model:  return updated object with stored parameter values (iFunc)
%         axes:   return the axes used for evaluation (cell of vector/matrix)
%         name:   return model name (char)
%
% ex:     b=feval(gauss,[1 2 3 4]); feval(gauss*lorz, [1 2 3 4, 5 6 7 8]);
%           feval(gauss,'guess', -5:5, -abs(-5:5))
%
% Version: $Date$
% See also iFunc, iFunc/fit, iFunc/plot

% handle input iFunc arrays
signal=[]; ax=[]; name='';

try
  inputname1=inputname(1);
  % handle bug in Matlab R2015-2016 for inputname shifted in feval
  % make sure argument is indeed what we expect
  if ~isa(evalin('caller', inputname1), 'iFunc') inputname1=''; end
catch
  inputname1 = '';
end


if isa(model, 'iFunc') && numel(model) > 1
  signal = {}; ax={}; name={};
  for index=1:numel(model)
    [signal{end+1}, model(index), ax{end+1}, name{end+1}] = feval(model(index), p, varargin{:});
  end
  if numel(signal) == 1, 
    signal=signal{1}; ax=ax{1};
  end
  signal = reshape(signal, size(model));
  ax     = reshape(ax,     size(model));
  name   = reshape(name,   size(model));
  return
end

% handle input parameter 'p' ===================================================

if nargin < 2
  p = [];
end

if ischar(model) && isa(p, 'iFunc')
  % call to an iFunc method
  signal = builtin('feval', model, p, varargin{:});
  return
end

if isa(p, 'iData')
  varargin = { p varargin{:} }; % will evaluate on iData axes with guesses pars
  p = NaN;
end

if iscell(p) && ~isempty(p) % as parameter cell (iterative function evaluation)
  signal = {}; ax={}; name={};
  for index=1:numel(p)
    [signal{end+1}, model, ax{end+1}, name{end+1}] = feval(model, p{index}, varargin{:});
  end
  if numel(signal) == 1, 
    signal=signal{1}; ax=ax{1}; name=name{1};
  end
  return
end

% some usual commands 
if strcmp(p, 'current'), p=model.ParameterValues; end
if ~isempty(p) && ischar(p)
  ax=[]; name=model.Name;
  if strcmp(p, 'plot')
    signal=plot(model);
    return
  elseif strcmp(p, 'identify')
    signal=evalc('disp(model)');
    return
  elseif ~strcmp(p, 'guess')
    p'
    disp([ mfilename ': Unknown parameter value in Model ' model.Name '. Using "guess" instead.'])
    p=[];
  end
elseif isa(p, 'iFunc')
  p=p.ParameterValues;
elseif ~isnumeric(p) && ~isempty(p) && ~isstruct(p)
  error([ 'iFunc:' mfilename ], [ 'Starting parameters "p" should be given as a vector, structure, character or empty, not ' class(p) ' length ' num2str(numel(p))]);
end

% convert a structure of parameters into an array matching model parameter
% names.
if isstruct(p)
  new = [];
  for index=1:length(model.Parameters)
    if isfield(p, model.Parameters{index})
      new = [ new p.(model.Parameters{index}) ];
    end
  end
  if length(new) == length(model.Parameters)
    p = new;
  else
    p
    disp([ 'Model ' model.Name ' parameters:' ])
    model.Parameters
    error([ 'iFunc:' mfilename ], 'Fields of the parameters "p" given as a structure do not match the model Parameters.');
  end
end

if ~ischar(p)
  p = p(:);
end

% handle varargin ==============================================================
% handle case where varargin contains itself model cell as 1st arg for axes and
% Signal

signal_in_varargin = []; % holds the index of a Signal after Axes in varargin
if ~isempty(varargin) 
  this = varargin{1};
  if iscell(this) % feval(model, p, {x y z ... Signal}, ...)
    Axes = this;
    if length(Axes) > abs(model.Dimension)
      Signal = Axes{abs(model.Dimension)+1};
      Axes   = Axes(1:abs(model.Dimension));
    end
    if ~isempty(Signal) && (ndims(Signal) == abs(model.Dimension) ...
      || (numel(Signal) == length(Signal) && numel(Signal) == numel(Axes{1})))
      Axes{end+1} = Signal; 
      signal_in_varargin = length(Axes);
    end
    varargin=[ Axes{:} varargin(2:end) ];
  elseif (isstruct(this) && isfield(this, 'Axes')) || isa(this, 'iData')
    % feval(model, p, struct('Axes','Signal'), ...)
    Signal = {};
    if isfield(this,'Signal')  
      Signal  = this.Signal;
      if isfield(this,'Monitor') && ~isa(this, 'iData')
        try
          Signal  = bsxfun(@rdivide,Signal, this.Monitor);
        end
      end
    end

    if isa(this, 'iData')
      Axes=cell(1,ndims(this));
      for index=1:ndims(this)
        Axes{index} = getaxis(this, index);
      end
    elseif isfield(this,'Axes')    Axes    = this.Axes; 
    end
    if ~isempty(Signal) && (ndims(Signal) == abs(model.Dimension) ...
      || (numel(Signal) == length(Signal) && numel(Signal) == numel(Axes{1})))
      Axes{end+1} = Signal; 
      signal_in_varargin = length(Axes);
    end
    varargin= [ Axes{:} varargin(2:end) ];
  end
  clear this Axes Signal
end

name=model.Name;
guessed = '';
% guess parameters ========================================================
[p,guessed,signal_in_varargin,ax] = iFunc_feval_guess_p(model, p, signal_in_varargin, varargin{:});

% format parameters and constraints as columns 
p = p(:);
if isfield(model.Constraint,'min')
  model.Constraint.min  =model.Constraint.min(:);
end
if isfield(model.Constraint,'max')
  model.Constraint.max  =model.Constraint.max(:);
end
if isfield(model.Constraint,'fixed')
  model.Constraint.fixed=model.Constraint.fixed(:);
end
if isfield(model.Constraint,'set')
  model.Constraint.set  =model.Constraint.set(:);
end

% apply constraints (fixed are handled in 'fits' -> forwarded to the optimizer)
if isfield(model.Constraint,'min')
  i = find(isfinite(model.Constraint.min));
  if ~isempty(i)
    p(i) = max(p(i), model.Constraint.min(i));
  end
end
if isfield(model.Constraint,'max')
  i = find(isfinite(model.Constraint.max));
  if ~isempty(i)
    p(i) = min(p(i), model.Constraint.max(i));
  end
end

% apply 'set' Constraints (with char)
p = iFunc_feval_set(model, p, varargin{:});

model.ParameterValues = p;

if ~isempty(inputname1)
  assignin('caller',inputname1,model); % update in original object
  if numel(ax) == model.Dimension+1
    signal = ax{end};
    ax(end) = [];
  end
elseif ~isempty(guessed)
  signal = model.ParameterValues;
end

% return here with syntax:
% feval(model) when model.ParameterValues is empty
% feval(model, 'guess')
if ~isempty(guessed)
  return
end

% guess axes ===================================================================
varargin = iFunc_feval_guess_axes(model, p, varargin{:});

% evaluate expression ==========================================================

% Eval contains both the Constraint and the Expression
% in case the evaluation is empty, we compute it (this should better have been done before)
if isempty(model.Eval) 
  model.Eval=cellstr(model);
  if ~isempty(inputname1)
    assignin('caller',inputname1,model); % update in original object
  end
end

% make sure we have enough parameter values wrt parameter names, else
% append 0's
if length(p) < length(model.Parameters)
  p = transpose([ p(:) ; zeros(length(model.Parameters) - length(p), 1) ]);
end

model.ParameterValues = p; % store current set of parameters

% request evaluation in sandbox, but should remove Signal after axes
if ~isempty(signal_in_varargin) && length(varargin) >= signal_in_varargin
  varargin(signal_in_varargin) = []; % remove Signal from arguments for evaluation (used in Guess)
  signal_in_varargin = [];
end
[signal,ax,p,model] = iFunc_feval_expr(model, varargin{:});

%model.ParameterValues = p; % store current set of parameters (updated)

p    = sprintf('%g ', p(:)'); if length(p) > 20, p=[ p(1:20) '...' ]; end
name = [ model.Name '(' p ') ' ];

if ~isempty(inputname1)
  assignin('caller',inputname1,model); % update in original object
end















% ==============================================================================
function [signal,iFunc_ax,p,this] = iFunc_feval_expr(this, varargin)
% private function to evaluate an expression in a reduced environment so that 
% internal function variables do not affect the result.

signal = [];
% assign parameters and axes for the evaluation of the expression, in case this is model char
% p already exists, we assign axes, re-assign varargin if needed
iFunc_ax = 'x y z t u v w ';
iFunc_dim = abs(this.Dimension);

if iFunc_dim && numel(varargin) >= iFunc_dim
  eval([ '[' iFunc_ax(1:(2*iFunc_dim)) ']=deal(varargin{' mat2str(1:iFunc_dim) '});' ]);
  % remove axes from varargin -> leaves additional optional arguments to the function
  varargin(1:iFunc_dim) = []; 
end



% EVALUATE now ...
% in the evaluation:
% * x,y,z,...        hold  the axes
% * p                holds the numerical values of the parameters (row)
% * struct_p         holds the parameters as a structure
% * this.Parameters holds the names of these parameters
% 

p       = reshape(this.ParameterValues,1,numel(this.ParameterValues));
% if we wish to have parameters usable as a structure
struct_p= cell2struct(num2cell(p),strtok(this.Parameters),2);

try
  this.Eval = cellstr(this.Eval);
  this.Eval = this.Eval(~strncmp('%', this.Eval, 1)); % remove comment lines
  eval(sprintf('%s\n', this.Eval{:}));
catch
  disp([ 'Error: Could not evaluate Expression in model ' this.Name ' ' this.Tag ]);
  disp(this)
  this.Eval
  lasterr
  save iFunc_feval_error
  error([ 'iFunc:' mfilename ], [ 'Failed model evaluation. Saved state in ' fullfile(pwd,'iFunc_feval_error') ]);
end

% copy the actual axes, in case they have been changed during evaluation
if nargout > 1 && iFunc_dim
  iFunc_ax = eval([ '{' iFunc_ax(1:(2*iFunc_dim)) '}' ]);
else
  iFunc_ax = [];
end

% ==============================================================================
function p = iFunc_feval_set(this, p, varargin)
% private function to evaluate a parameter set expression in a reduced environment so that 
% internal function variables do not affect the result.
  if isfield(this.Constraint, 'set')
    i = find(~cellfun('isempty', this.Constraint.set)); i=i(:)';
  else i=[]; end
  if ~isempty(i)

    ax = 'x y z t u v w';
    this.Dimension = abs(this.Dimension);
    
    if this.Dimension
      eval([ '[' ax(1:(2*this.Dimension)) ']=deal(varargin{' mat2str(1:this.Dimension) '});' ]);
    end
    if length(varargin) > this.Dimension && ~isempty(varargin{this.Dimension+1}) && isnumeric(varargin{this.Dimension+1})
      signal = varargin{this.Dimension+1};
    else
      signal = 1;
    end
    clear ax

    for index=i
      try
        if isa(this.Constraint.set{index}, 'function_handle') && ...
           nargout(this.Constraint.set{index}) == 1
          n = nargin(this.Constraint.set{index});
          if n > 0 && length(varargin) >= n
            p(index) = feval(this.Constraint.set{index}, p, varargin{1:n});
          else
            p(index) = feval(this.Constraint.set{index}, p, varargin);
          end
        elseif ischar(this.Constraint.set{index})
          p(index) = eval(this.Constraint.set{index});
        end
      catch ME
        disp([ mfilename ': Constraints: ' ME.message ])
        warning([ 'iFunc:' mfilename ], 'Could not evaluate model.Constraint.set on p(%i):', index);
        disp(this.Constraint.set{index})
      end % try
      
    end
  end

