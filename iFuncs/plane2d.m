function signal=plane2d(p, x, y, signal)
% signal = plane2d(p, x, y, {signal}) : Planar function
%
%   iFunc/plane2d Planar function (fit 2D function/model)
%       signal = p(1)*x+p(2)*y+p(3)
%   the function called with a char argument performs specific actions.
%
% input:  p: plane2d model parameters (double array)
%            p = [  'Slope_X' 'Slope_Y' 'Background' ]
%          or action e.g. 'identify', 'guess', 'plot' (char)
%         x: axis along signal rows    (double)
%         y: axis along signal columns (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value or information structure ('guess', 'identify','plot')
% ex:     signal=plane2d([1 2 .5 .2 .3 30 .2], -2:.1:2, -3:.1:3); or p=plane2d('guess',x,y,signal);
%
% Version: $Revision: 1.1 $
% Reference: http://en.wikipedia.org/wiki/Gaussian_function
% See also iData, ifitmakefunc, gauss, iData/fits

% 2D function created by ifitmakefunc from template2d (iFit/iFuncs)
%   Please retain the function definition structure as defined below
%   in most cases, just fill-in the information when HERE is indicated

  if nargin >= sum(isstrprop('x, y','alpha'))+2 && isnumeric(signal) && ~isempty(signal) ...
                     && isnumeric(x) && ~isempty(x) ...
                     && isnumeric(y) && ~isempty(y)
  %   guess: model('guess', x, y,signal)
  %   guess: model(p,       x, y,signal)
    signal = guess(x, y,signal);
  elseif nargin >= sum(isstrprop('x, y','alpha'))+1 && isnumeric(p) && ~isempty(p) ...
                 && isnumeric(x) && ~isempty(x) ...
                 && isnumeric(y) && ~isempty(y)
  %   evaluate: model(p,x, y, ...)
    signal = evaluate(p, x, y);
  elseif nargin >= sum(isstrprop('x, y','alpha'))+1 && isnumeric(p) && isnumeric(x) && isnumeric(y) ...
                     && any(size(p) == size(y)) && any(size(x) == size(y))
  %   guess: model(x, y,signal) with numel(x, y)==size(signal)
    signal = guess(p,x, y);
  elseif (nargin >= 3 && isnumeric(p) && ~isempty(p) && (isempty(x) || isempty(y))) ...
      || (nargin >= 2 && isnumeric(p) && ~isempty(p) && isempty(x))
  %   evaluate: model(p,[])
    signal = feval(mfilename, p);
  elseif nargin == sum(isstrprop('x, y','alpha'))+1 && isempty(p) && isnumeric(x) && ~isempty(x) && ~isempty(y)
  %   identify: model([],x, y)
    info = identify; signal=[];
    info.Axes   = { x, y };
    info.Guess  = iFuncs_private_guess(info.Axes, signal, info.Parameters);
    info.Values = evaluate(info.Guess, info.Axes{:});
    signal = info;
  elseif nargin == 1 && isnumeric(p) && ~isempty(p) 
  %   identify: model(p)
    signal = identify;
    signal.Guess  = p;
    width1    = 0; width2    = 0; 
    position1 = 0; position2 = 0;
    parameter_names = lower(signal.Parameters);
    % HERE default axes to represent the model when parameters are given <<<<<<<
    % test parameter names
    for index_p=1:length(parameter_names)
      if  ~isempty(strfind(parameter_names{index_p}, 'width')) ...
        | ~isempty(strfind(parameter_names{index_p}, 'tau')) ...
        | ~isempty(strfind(parameter_names{index_p}, 'damping'))
        width = width+abs(p(index_p));                     % kind of embrace all peaks
      elseif ~isempty(strfind(parameter_names{index_p}, 'centre')) ...
        |    ~isempty(strfind(parameter_names{index_p}, 'center')) ...
        |    ~isempty(strfind(parameter_names{index_p}, 'position'))
        if position == 0, position = p(index_p);
        else position = mean( [position p(index_p)] ); end % kind of center on all peaks
      end
    end
    if width == 0
      if position == 0, width=1; else width=abs(position/2); end
    end
    % we use a single axis, repeated along dimensions as it is not possible to identify
    % how to set axes from only the parameter names.
    x = linspace(position1-3*width1,position1+3*width1, 50);
    signal.Axes = cell(1, sum(isstrprop('x, y','alpha')));
    for index=1:sum(isstrprop('x, y','alpha'))
      signal.Axes{index}   =  x;
    end
    signal.Values = evaluate(signal.Guess, signal.Axes{:});
  elseif nargin == 1 && ischar(p) && strcmp(p, 'plot')
    signal = feval(mfilename);
    if signal.Dimension == 2
      h=surf(signal.Axes{2}, signal.Axes{1}, signal.Values); set(h,'Edgecolor','none');
      xlabel('columns (2nd axis)'); ylabel('rows (1st axis)'); % plot convention reversed wrt matrix
    elseif signal.Dimension == 3
      isosurface(signal.Axes{2}, signal.Axes{1}, signal.Axes{3}, signal.Values, median(signal.Values(:)));
    else
      return
    end
    title(mfilename);
  elseif nargin == 0
    x = linspace(-2,2, 100);
    signal.Axes = cell(1, sum(isstrprop('x, y','alpha')));
    for index=1:sum(isstrprop('x, y','alpha'))
      signal.Axes{index}   =  x;
    end
    signal = feval(mfilename, [], signal.Axes{:});
  else
    signal = identify;
  end

end
%                                        end of model main gauss2d
% ------------------------------------------------------------------------------

% inline: evaluate: compute the model values
function signal = evaluate(p, x, y)
  if isempty(x) || isempty(y) || isempty(p), signal=[]; return; end
  % check if x,y are not mesh-like, then have to be computed on a mesh
  if isvector(x) && isvector(y)
    x=x(:); y=y(:)';
    [x, y]=ndgrid(x, y);
  end
  
  % HERE is the model evaluation <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    signal = p(1)*x+p(2)*y+p(3);

  signal = reshape(signal, size(x));
end

% inline: identify: return a structure which identifies the model
function signal =identify()
  % HERE are the parameter names
  parameter_names = { 'Slope_X' 'Slope_Y' 'Background' };
  %
  signal.Type           = 'iFit fitting function';
  signal.Name           = [ 'Planar function (2D) [' mfilename ']' ];
  signal.Parameters     = parameter_names;
  signal.Dimension      = 2;         % dimensionality of input space (axes) and result
  signal.Guess          = [];        % default parameters
  signal.Axes           = {};        % the axes used to get the values
  signal.Values         = [];        % default model values=f(p)
  signal.function       = mfilename;
end

% inline: guess: guess some starting parameter values and return a structure
function info=guess(x, y,signal)
  info       = identify;  % create identification structure
  info.Axes  = { x, y };
  % fill guessed information
  info.Guess = iFuncs_private_guess(info.Axes, signal, info.Parameters);
  info.Values= evaluate(info.Guess, info.Axes{:});
end


