function [b, info] = ieval(a, model, pars, varargin)
% [b, info] = ieval(a, model, varargin) evaluate a function on the axes of an object
%
%   @iData/ieval applies the function 'model' using the axes of the object 'a'
%     and function parameters 'pars' with additional parameters.
%     The model function must follow the syntax:
%       y = model(pars, axis1, axis2, ...)
%     If the model function is specified as a cellstr containing references to
%     lower dimensionality functions as the iData object, the model is build
%     as the product of these functions. This enables to describe a N-D data set
%     as series of e.g. 1D functions.
%
% input:  a: object or array (iData)
%         model: model function (char/function handle/cellstr)
%         pars: model parameters (double array)
%         additional parameters may be passed.
% output: b: result of evaluation (iData)
%         info: structure giving information about the model
% ex:     b=ieval(a,'gauss',[1 2 3 4]); or ieval(a, {'gauss','lorentz'}, [1 2 3 4, 5 6 7 8]);
%
% See also iData, feval
% Contributed code (Matlab Central): 
%   genop: Douglas M. Schwarz, 13 March 2006

% private functions: 
%   genop: Douglas M. Schwarz, 13 March 2006

if nargin < 2
  iData_private_error(mfilename, 'Function evaluation requires at least 3 parameters');
end
if nargin < 3
  pars = [];
end
if nargin < 4
  varargin = {};
end

if length(a) > 1
  b = a;
  for index=1:length(a(:))
    b(index) = ieval(a(index), model, pars, varargin);
  end
  b = reshape(b, size(a));
  return
end

Axes   = cell(1,ndims(a));
for index=1:ndims(a)
  Axes{index} = getaxis(a, index);  % loads object axes, or 1:end if not defined 
end

% evaluate model signal
if ischar(model) | isa(model, 'function_handle')
  if nargout > 1 | isempty(pars)
    info = feval(model,'identify');  % get identification info
  end
  if isempty(pars), pars=info.Guess; end
  if ~isempty(varargin)
    Model = feval(model, pars, Axes{:}, varargin{:});
  else
    Model = feval(model, pars, Axes{:});
  end 
elseif iscell(model)
  % identify the model dimensionality
  axis_index=1;
  pars_index=1;
  model_values={};
  model_ndims ={};
  model_pars=[]; model_namepars={}; model_name='';
  for index=1:length(model(:))
    model_info = feval(model{index},'identify');  % get identification info
    % check dimensions: are enough axes and parameters available ?
    if length(Axes(:)) < axis_index+model_info.Dimension-1
      iData_private_error([ mfilename '/' model{index} ], ...
        [ 'Axis length is ' num2str(length(Axes(:))) ' but the axis ' ...
          num2str(axis_index+model_info.Dimension-1) ' is requested in ' ...
          num2str(index) '-th model function ' char(model{index}) ' when fitting object ' a.Tag ]);
    end
    if ~isempty(pars) & length(pars) < pars_index+length(model_info.Parameters)-1
      iData_private_error([ mfilename '/' model{index} ], ...
        [ 'Parameters length is ' num2str(length(pars)) ' but the parameter ' ...
          num2str(pars_index+length(model_info.Parameters)-1) ' is requested in ' ...
          num2str(index) '-th model function ' char(model{index}) ' when fitting object ' a.Tag ]);
    end
    if isempty(pars)
      this_pars = model_info.Guess;
    else
      this_pars = pars(pars_index:(pars_index+length(model_info.Parameters)-1));
    end
    this_axes = Axes{axis_index:(axis_index+model_info.Dimension-1)};
    model_pars = [ model_pars this_pars ];
    % evaluate sub-model
    if ~isempty(varargin)
      model_value  = feval(model{index}, this_pars, this_axes, varargin{:} );
    else
      model_value  = feval(model{index}, this_pars, this_axes);
    end
    model_value  = squeeze(model_value);
    model_values = { model_values{:} ;  model_value }; % append to model values
    model_name   = [ model_name '*' func2str(model{index}) ];
    
    % get the dimensionality of sub-model
    n = size(model_value);
    if     all(n == 0), n=0; continue;
    else n = n(find(n > 1)); end
  
    % assign individual model dimensions
    model_ndim  =ones(1,ndims(a));
    model_ndim(axis_index:(axis_index+model_info.Dimension-1)) = n;
    model_ndims ={ model_ndims{:} ; model_ndim };
    axis_index  =axis_index+model_info.Dimension;
    pars_index  =pars_index+length(model_info.Parameters);
    model_namepars  = { model_namepars{:} ; model_info.Parameters };
  end % for sub-models
  % now make up the product of sub-space models
  Model = model_values{1};
  for index=2:length(model_values)
    Model = genop(@times, Model, model_values{index});
  end
  
  info.Type      = 'iFit fitting function';
  info.Name      = model_name;
  info.Parameters= model_namepars;
  info.Dimension = sum(cell2mat(model_ndims));
  info.Guess     = model_pars;
else
  iData_private_error(mfilename, ['2nd argument "model" must be specified as a char, a cellstr or function handles. Currently of type ' class(model) ]);
end

% build the output iData object
b = copyobj(a);
if isnumeric(Model)
  setalias(b,'Signal', Model);
  setalias(b,'Error', 0);
end
b = iData_private_history(b, mfilename, a, model, pars, varargin{:});  
% final check
b = iData(b);




