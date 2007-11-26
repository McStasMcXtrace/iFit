function b = feval(a, model, pars, varargin)
% b = feval(a, model, varargin) evaluate a function on the axes of an object
%
%   @iData/feval applies the function 'model' using the axes of the object 'a'
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
% ex:     b=feval(a,'gauss',[1 2 3 4]); or feval(a, {'gauss','lorentz'}, [1 2 3 4, 5 6 7 8]);
%
% See also iData, feval
% Contributed code (Matlab Central): 
%   genop: Douglas M. Schwarz, 13 March 2006

% private functions: 
%   genop: Douglas M. Schwarz, 13 March 2006

if nargin < 3
  iData_private_error(mfilename, 'Function evaluation requires at least 3 parameters');
end

Axes   = cell(1,ndims(a));
for index=1:ndims(a)
  Axes{index} = getaxis(a, index);  % loads object axes, or 1:end if not defined 
end

% evaluate model signal
if ischar(model) | isa(model, 'function_handle')
  Model = feval(model, pars, Axes{:}, varargin{:});
elseif iscell(model)
  % identify the model dimensionality
  axis_index=1;
  pars_index=1;
  model_values={};
  model_ndims ={};
  for index=1:length(model(:))
    model_info = feval(model{index},'identify');  % get identification info
    % check dimensions: are enought axes and parameters available ?
    if length(a_axes(:)) < axis_index+model_info.Dimension-1
      iData_private_error([ mfilename '/' model{index} ], ...
        [ 'Axis length is ' num2str(length(a_axes(:))) ' but the axis ' ...
          num2str(axis_index+model_info.Dimension-1) ' is requested in ' ...
          num2str(index) '-th model function ' char(model{index}) ' when fitting object ' a.Tag ]);
    elseif length(pars) < pars_index+length(model_info.Parameters)-1
      iData_private_error([ mfilename '/' model{index} ], ...
        [ 'Parameters length is ' num2str(length(pars)) ' but the parameter ' ...
          num2str(pars_index+length(model_info.Parameters)-1) ' is requested in ' ...
          num2str(index) '-th model function ' char(model{index}) ' when fitting object ' a.Tag ]);
    end
    % evaluate sub-model
    model_value  = feval(model{index}, ...
      pars(pars_index:(pars_index+length(model_info.Parameters)-1)), ...
      a_axes{axis_index:(axis_index+model_info.Dimension-1)}, varargin{:} );
    model_value  = squeeze(model_value);
    model_values = { model_values{:} ;  model_value }; % append to model values
    
    % get the dimensionality of sub-model
    n = size(model_value);
    if     all(n == 0), n=0;
    else n = n(find(n > 1)); end
  
    % assign individual model dimensions
    model_ndim  =ones(1,ndims(a));
    model_ndim(axis_index:(axis_index+model_info.Dimension-1)) = n;
    model_ndims ={ model_ndims{:} ; model_ndim };
    axis_index  =axis_index+model_info.Dimension;
    pars_index  =pars_index+length(model_info.Parameters);
  end % for sub-models
  % now make up the product of sub-space models
  Model = model_values{1};
  for index=2:length(model_values)
    Model = genop(@times, Model, model_values{index});
  end
else
  iData_private_error(mfilename, ['2nd argument "model" must be specified as a char, a cellstr or function handles. Currently of type ' class(model) ]);
end

% build the output iData object
b = copyobj(a);
setalias(b,'Signal', Model);
b = iData_private_history(b, mfilename, a, model, pars, varargin{:});  
% final check
b = iData(b);



