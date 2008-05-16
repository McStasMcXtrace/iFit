function [b, Info] = ieval(a, model, pars, varargin)
% [b, Info] = ieval(a, model, varargin) evaluate a function on the axes of an object
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
%         Info: structure giving information about the model
% ex:     b=ieval(a,'gauss',[1 2 3 4]); or ieval(a, {'gauss','lorentz'}, [1 2 3 4, 5 6 7 8]);
%
% Contributed code (Matlab Central): 
%   genop: Douglas M. Schwarz, 13 March 2006
%
% See also iData, feval

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
  Info = feval(model,'identify');  % get identification Info
  % check dimensionality
  if isfield(Info, 'Dimension')
    if Info.Dimension ~= ndims(a)
      if rem(ndims(a), Info.Dimension) == 0
        list={};
        for index=1:(ndims(a)/Info.Dimension)
          list = { list{:} ; model };
        end
        [b, Info] = ieval(a, list, pars);
      else
        iData_private_error([ mfilename '/' model ], ...
          [ 'Can not build ' num2str(ndims(a)) ' dimensionality model from ' ...
            'model function ' char(model) ' of dimension ' num2str(Info.Dimension) ' when fitting object ' a.Tag ]);
      end
      % check Info fields
      if ~isfield(Info,'Type') Info.Type = 'iFit fitting function'; end
      if ~isfield(Info,'Name') Info.Name = char(model); end
      if ~isfield(Info,'Parameters') & length(pars)
        Info.Parameters = {};
        for index=1:length(pars), Info.Parameters{index} =['Parameter_' num2str(index)]; end
      end
      if ~isfield(Info,'Guess') Info.Guess= rand(1,length(pars)); end % default parameters
      return
    end
  end
  if isfield(Info, 'Guess') & isempty(pars), pars=Info.Guess; end
  if ~isempty(varargin)
    Model = feval(model, pars, Axes{:}, varargin{:});
  elseif ~isempty(pars)
    Model = feval(model, pars, Axes{:});
  else
    try
    Model = feval(model, Axes{:});
    catch
    Model = feval(model, double(a));
    end
  end 
elseif iscell(model)
  % identify the model dimensionality
  axis_index=1;
  pars_index=1;
  model_values={};
  model_ndims ={};
  model_pars=[]; model_namepars={}; model_name='';
  for index=1:length(model(:))
    model_Info = feval(model{index},'identify');  % get identification Info
    % check Info fields
    if ~isfield(model_Info,'Type') model_Info.Type = 'iFit fitting function'; end
    if ~isfield(model_Info,'Name') model_Info.Name = char(model); end
    if ~isfield(model_Info,'Guess') model_Info.Guess= rand(1,length(model_Info.Parameters)); end 
    % check dimensions: are enough axes and parameters available ?
    if length(Axes(:)) < axis_index+model_Info.Dimension-1
      iData_private_error([ mfilename '/' model{index} ], ...
        [ 'Axis length is ' num2str(length(Axes(:))) ' but the axis ' ...
          num2str(axis_index+model_Info.Dimension-1) ' is requested in ' ...
          num2str(index) '-th model function ' char(model{index}) ' when fitting object ' a.Tag ]);
    end
    if isnumeric(pars) & ~isempty(pars) & length(pars) < pars_index+length(model_Info.Parameters)-1
      iData_private_error([ mfilename '/' model{index} ], ...
        [ 'Parameters length is ' num2str(length(pars)) ' but the parameter ' ...
          num2str(pars_index+length(model_Info.Parameters)-1) ' is requested in ' ...
          num2str(index) '-th model function ' char(model{index}) ' when fitting object ' a.Tag ]);
    end
    if isempty(pars) | strcmp(pars,'identify')
      this_pars = model_Info.Guess;
    else
      this_pars = pars(pars_index:(pars_index+length(model_Info.Parameters)-1));
    end
    this_axes = Axes{axis_index:(axis_index+model_Info.Dimension-1)};
    model_pars = [ model_pars this_pars ];
    % evaluate sub-model
    if ~isempty(varargin)
      model_value  = feval(model{index}, this_pars, this_axes, varargin{:} );
    else
      model_value  = feval(model{index}, this_pars, this_axes);
    end
    model_value  = squeeze(model_value);
    model_values = { model_values{:} ;  model_value }; % append to model values
    fname = [ char(model{index}) '_f' num2str(index) ];
    if isempty(model_name)
      model_name   = fname;
    else
      model_name   = [ model_name '*' fname ];
    end
    
    % get the dimensionality of sub-model
    n = size(model_value);
    if     all(n == 0), n=0; continue;
    else n = n(find(n > 1)); end
  
    % assign individual model dimensions
    model_ndim  =ones(1,ndims(a));
    model_ndim(axis_index:(axis_index+model_Info.Dimension-1)) = n;
    model_ndims ={ model_ndims{:} ; model_ndim };
    axis_index  =axis_index+model_Info.Dimension;
    pars_index  =pars_index+length(model_Info.Parameters);
    pnames = model_Info.Parameters;
    pnames(2,:) = { ['_f' num2str(index) ] };
    pnames      = strcat(pnames(1,:), pnames(2,:));  
    model_namepars  = { model_namepars{:} , pnames{:} };
  end % for sub-models
  % now make up the product of sub-space models
  if ~isempty(model_values)
    Model = model_values{1};
    for index=2:length(model_values)
      Model = genop(@times, Model, model_values{index});
    end
  else Model=[]; model_pars=Info.Guess;
  end
  
  Info.Type      = 'iFit fitting function';
  Info.Name      = model_name;
  Info.Parameters= model_namepars;
  Info.Dimension = sum(cell2mat(model_ndims));
  Info.Guess     = model_pars;
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

