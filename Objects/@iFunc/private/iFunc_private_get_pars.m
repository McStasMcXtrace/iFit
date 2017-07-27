function [pars, pars_isstruct, ax] = iFunc_private_get_pars(model, pars, varargin)
  % build a parameter vector from starture, char or array
  % fills all non specified parameters
  %
  % pars: the parameter vector with all the required values
  % pars_isstruct: vector indicating which parameters where given as a structure

  if nargin < 3, ax = {}; 
  else ax = varargin; end

  pars_isstruct=[];
  if isempty(pars), pars='current'; end
  if ischar(pars) && strcmp(pars,'current')
    pars = model.ParameterValues;
  end
  if ischar(pars) && ~strcmp(pars,'guess')
    pars = str2struct(pars);
  end
  
  p0 = model.ParameterValues;
  % test for various parameter input
  if (isnumeric(pars) && length(pars) < length(model.Parameters)) ...
    || isstruct(pars)
    if isempty(p0) % current values are unset
      [p0, model, ax] = feval(model, 'current', varargin{:});
      p0 = model.ParameterValues;
    end
  end
  
  if strcmp(pars,'guess')
    % force guess
    [p0, model, ax] = feval(model, 'guess', varargin{:});
    pars = model.ParameterValues;
  elseif isnumeric(pars) && length(pars) < length(model.Parameters)
    % incomplete parameter array
    pars = [ pars(:)' p0((numel(pars)+1):end)' ];
    
  elseif isstruct(pars) % structure with possibly not all parameters specified
    if isempty(fieldnames(pars))
      pars_isstruct = 1:numel(model.Parameters);
      pars = p0;
    else
      % we scan the given structure members. When matches a name, we set the p0 value
      f_given=fieldnames(pars);
      f_model=model.Parameters;
      for index=1:length(f_model)  
        match = strcmp(strtok(f_model{index}), f_given);
        if ~isempty(match) && match
          pars_isstruct = [ pars_isstruct index ];
          p0(index) = pars.(f_given{match});
        end
      end
      pars = p0;
    end
  end
  pars = reshape(pars, [ 1 numel(pars)]); % a single row
  
  % get axes
  if nargout>2 && isempty(ax) || numel(ax) < model.Dimension
    [~, ~, ax] = feval(model, pars, ax{:});
  end
  if numel(ax) > model.Dimension
    ax = ax(1:model.Dimension);
  end

