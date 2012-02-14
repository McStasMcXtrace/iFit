function [b, Info] = ieval(a, model, pars, varargin)
% [b, Info] = ieval(a, model, pars, ...) evaluate a function on the axes of an object
%
%   @iData/ieval applies the function 'model' using the axes of the object 'a'
%     and function parameters 'pars' with optional additional parameters.
%     The model function must follow the syntax:
%         y = model(pars, axis1, axis2, ...)
%     If the model function is specified as a cellstr containing references to
%       lower dimensionality functions as the iData object, the model is built
%       as the product of these functions. This enables to describe a N-D data set
%       as series of e.g. 1D functions.
%     The special call to ieval(a,model,'guess') and ieval(a,model,'identify') will 
%       request a model parameter guess and identification resp.
%
% input:  a: object or array (iData)
%         model: model function (char/function handle/cellstr)
%         pars: model parameters (vector) or 'guess' or 'identify'
%         additional parameters may be passed.
% output: b: result of evaluation (iData)
%         Info: structure giving information about the model
% ex:     b=ieval(a,'gauss',[1 2 3 4]); ieval(a, {'gauss','lorentz'}, [1 2 3 4, 5 6 7 8]);
%           ieval(a,'gauss','guess')
%
% Version: $Revision: 1.20 $
% See also iData, feval, iData/fits

% private functions: 
%   genop: Douglas M. Schwarz, 13 March 2006

if nargin < 1
  iData_private_error(mfilename, 'Function evaluation requires at least 1 parameter (iData)');
end
if nargin < 2
  model = 'gauss';
end
if nargin < 3
  pars = [];
end

% handle array of objects to fit iteratively
if length(a) > 1
  b = a;
  for index=1:length(a(:))
    b(index) = ieval(a(index), model, pars, varargin{:});
  end
  b = reshape(b, size(a));
  return
end

% get the object axes
Axes   = cell(1,ndims(a));
for index=1:ndims(a)
  Axes{index} = getaxis(a, index);  % loads object axes, or 1:end if not defined 
end

% evaluate model signal
if ischar(model) | isa(model, 'function_handle') % a single model ==============

  Info = feval(model,'identify');  % get identification Info, esp. for Dimension
  
  % check dimensionality
  if isfield(Info, 'Dimension')
    if Info.Dimension <= ndims(a) % model dimensionality does not match object one
      if rem(ndims(a), Info.Dimension) == 0 % must be able to fill object dimensions with model
        list={};
        for index=1:(ndims(a)/Info.Dimension)
          list = { list{:} ; model };
        end
        model = list;
      else
        iData_private_error([ mfilename '/' char(model) ], ...
          [ 'Can not build ' num2str(ndims(a)) ' dimensionality model from ' ...
            'model function ' char(model) ' of dimension ' num2str(Info.Dimension) ...
            ' when fitting object ' a.Tag ]);
      end
    end
  end
  % and will call the model(cell) syntax below
end
   
if iscell(model) % a set of models (product) ===================================
  % identify the model dimensionality
  axis_index=1;
  pars_index=1;
  model_values={};
  model_ndims ={};
  model_pars=[]; model_namepars={}; model_name='';
  for index=1:length(model(:))
    % get identification Info
    model_Info = feval(model{index}, 'identify');

    % check dimensions: are enough axes and parameters available ?
    if length(Axes(:)) < axis_index+model_Info.Dimension-1
      iData_private_error([ mfilename '/' char(model{index}) ], ...
        [ 'Axis length is ' num2str(length(Axes(:))) ' but the axis ' ...
          num2str(axis_index+model_Info.Dimension-1) ' is requested in ' ...
          num2str(index) '-th model function ' char(model{index}) ' when fitting object ' a.Tag ]);
    end
    if isnumeric(pars) & ~isempty(pars) & length(pars) < pars_index+length(model_Info.Parameters)-1
      iData_private_error([ mfilename '/' char(model{index}) ], ...
        [ 'Parameters length is ' num2str(length(pars)) ' but the parameter ' ...
          num2str(pars_index+length(model_Info.Parameters)-1) ' is requested in ' ...
          num2str(index) '-th model function ' char(model{index}) ' when fitting object ' a.Tag ]);
    end
    
    % Axes for this model
    this_axes  = Axes(axis_index:(axis_index+model_Info.Dimension-1));
    
    % Signal for this model if guess needed
    if isempty(pars) || (ischar(pars) && any(strcmp(pars,{'guess','identify'})))
      Signal = get(a,'Signal'); Signal(~isfinite(Signal)) = 0;
      m = get(a,'Monitor'); m=real(m(:));
      if not(all(m == 1 | m == 0)),
        Signal = genop(@rdivide,Signal,m);
      end
      % determine indexes on which sum is required. The final Model
      % must match model_Info.Dimension
      j = ones(1,length(Axes(:)));
      j(axis_index+model_Info.Dimension-1) = 0; % these will not be summed-up
      if ~isvector(Signal)
        for i=length(j):-1:1
          if j(i), Signal=trapz(Signal, i)/size(Signal,i); end  % reduce Model dimensionality
        end
      end
    else Signal=[]; end
    
    if nargin < 4
      varargin = { Signal };
    elseif ~isempty(Signal)
      varargin = { Signal, varargin{:} };
    end 
    clear Signal
    
    if isempty(pars) || (ischar(pars) && any(strcmp(pars,{'guess','identify'})))
      % obtained guessed parameters if pars is missing or explicitely requested
      this_pars = feval(model{index}, 'guess', this_axes{:}, varargin{:});
      this_pars = this_pars.Guess;
    else
      this_pars = pars(pars_index:(pars_index+length(model_Info.Parameters)-1));
    end
    
    model_pars = [ model_pars this_pars ];
    
    if ~isempty(varargin)
      model_value  = feval(model{index}, this_pars, this_axes{:}, varargin{:} );
    else
      model_value  = feval(model{index}, this_pars, this_axes{:});
    end
    model_value  = squeeze(model_value);
    model_values = { model_values{:} ;  model_value }; % append to model values
    fname = char(model{index});
    if length(model(:)) > 1
      fname = [ fname '_f' num2str(index) ];
    end
    if isempty(model_name)
      model_name   = fname;
    else
      model_name   = [ model_name '*' fname ];
    end
    
    if isnumeric(model_value)
      % get the dimensionality of sub-model
      n = size(model_value);
      if   all(n == 0), n=0; continue;
      else n = n(find(n > 1)); end
    
      % assign individual model dimensions
      model_ndim  =ones(1,ndims(a));
      model_ndim(axis_index:(axis_index+model_Info.Dimension-1)) = n;
    elseif isstruct(model_value) && ischar(pars) && any(strcmp(pars,{'guess','identify'}))
      model_ndim   = model_value.Dimension;
      moldel_value = model_value.Values;
    end
    clear model_value
    
    model_ndims ={ model_ndims{:} ; model_ndim };
    axis_index  =axis_index+model_Info.Dimension;
    pars_index  =pars_index+length(model_Info.Parameters);
    pnames = model_Info.Parameters;
    if length(model(:)) > 1
      pnames(2,:) = { ['_f' num2str(index) ] };
      pnames      = strcat(pnames(1,:), pnames(2,:));
    end
    model_namepars  = { model_namepars{:} , pnames{:} };
  end % for sub-models
  
  % now make up the product of sub-space models
  if ~isempty(model_values)
    Model = model_values{1};
    for index=2:length(model_values)
      % sub-spaces above 1 must be normalized, not to affect the total intensity
      this_value = model_values{index};
      this_value = this_value/norm(this_value(:));
      Model = genop(@times, Model, this_value);
    end
  else Model=[]; model_pars=this_pars;
  end

  Info.Type      = 'iFit fitting function';
  Info.Name      = model_name;
  Info.Parameters= model_namepars;
  Info.Dimension = sum(cell2mat(model_ndims));
  Info.Guess     = model_pars;
  pars = model_pars;
else
  iData_private_error(mfilename, ['2nd argument "model" must be specified as a char, a cellstr or function handles. Currently of type ' class(model) ]);
end

% now got [Model, Info] for single model or cell

% assign return values
if isstruct(Model)
  Info = Model;
  pars  = Info.Guess;
  Model = Info.Values;
  if isempty(Model), 
    b=[]; return; 
  end
end

% build the output iData object
cmd=a.Command;
b = copyobj(a);
if isnumeric(Model)
  if iscellstr(model)
    model = char(strcat(strcat(model(1:(end-1)),'*'), model(end)))';
    model = model(:)';
  end
  m = get(b,'Monitor'); m=real(m); m=m(:);
  if not(all(m == 1 | m == 0)) && (numel(m) == 1 || numel(m) == numel(Model)),
    Model = Model.*m;
  end
  clear m
  setalias(b,'Signal', Model, char(Info.Name));
  b.Title = [ char(Info.Name) '(' b.Title ')' ];
  b.Label = b.Title;
  setalias(b,'Error', 0);
  setalias(b,'Parameters', pars, [ char(Info.Name) ' model parameters for ' char(a) ]);
  Info.Guess = pars;
  b.Data.Model = Info;
  setalias(b,'Model','Data.Model',[ char(Info.Name) ' model description for ' char(a) ]);
end
b.Command=cmd;
b = iData_private_history(b, mfilename, a, char(Info.Name), pars, varargin{:});
% final check
b = iData(b);

% end of ieval
