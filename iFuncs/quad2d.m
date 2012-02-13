function signal=quad2d(p, x, y, signal)
% signal = quad2d(p, x, y, {signal}) : 2D Quadratic function
%
%   iFunc/quad2d 2D Quadratic function (fit 2D function/model)
%     x0=p(2); y0=p(3); sx=p(4); sy=p(5); theta=p(6) [given in deg]
%     a = cos(theta)^2/2/sx/sx + sin(theta)^2/2/sy/sy;
%     b =-sin(2*theta)/4/sx/sx + sin(2*theta)/4/sy/sy;
%     c = sin(theta)^2/2/sx/sx + cos(theta)^2/2/sy/sy;
%     signal = p(1)*(a*(x-x0).^2+2*b*(x-x0).*(y-y0)+c*(y-y0).^2)+ p(7);
%   the function called with a char argument performs specific actions.
%
% input:  p: quad2d model parameters (double)
%            p = [  'Amplitude' 'Centre_X' 'Center_Y'
%                   'HalfWidth_X' 'HalfWidth_Y' 'Angle' 'Background' ] as a numerical array
%            the rotation angle is given in degrees.
%          or action e.g. 'identify', 'guess', 'plot' (char)
%         x: axis along rows    (double)
%         y: axis along columns (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value or information structure ('guess', 'identify','plot')
% ex:     signal=quad2d([1 2 .5 .2 .3 30 .2], -2:.1:2, -3:.1:3); or p=quad2d('guess',x,y,signal);
%
% Version: $Revision: 1.2 $
% Reference: http://en.wikipedia.org/wiki/Quadratic_function
% See also iData, ifitmakefunc, gauss, iData/fits

  if nargin >= 3 && isnumeric(p) && ~isempty(p) ...
                 && isnumeric(x) && ~isempty(x) ...
                 && isnumeric(y) && ~isempty(y)
  %   evaluate: model(p,x, y, ...)
    signal = evaluate(p, x, y);
  elseif nargin == 4 && isnumeric(signal) && ~isempty(signal) ...
                     && isnumeric(x) && ~isempty(x) ...
                     && isnumeric(y) && ~isempty(y)
  %   guess: model('guess', x,y,signal)
  %   guess: model(p,       x,y,signal)
    signal = guess(x,y,signal);
  elseif nargin == 3 && isnumeric(p) && isnumeric(x) && isnumeric(y) ...
                     && any(size(p) == size(y)) && any(size(x) == size(y))
  %   guess: model(x,y,signal) with numel(x,y)==size(signal)
    signal = guess(p,x,y);
  elseif (nargin == 3 && isnumeric(p) && ~isempty(p) && (isempty(x) || isempty(y))) ...
      || (nargin == 2 && isnumeric(p) && ~isempty(p) && isempty(x))
  %   evaluate: model(p,[])
    signal = feval(mfilename, p);
  elseif nargin == 3 && isempty(p) && isnumeric(x) && ~isempty(x) && ~isempty(y)
  %   identify: model([],x,y)
    info = identify; 
    info.Guess  = [ 1 mean(x(:)) mean(y(:)) std(x(:)) std(y(:)) 30 0 ];
    info.Axes   = { x y };
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
        if (~width1) width1 = width1+abs(p(index_p));                     % kind of embrace all peaks
        else         width2 = width2+abs(p(index_p)); end
      elseif ~isempty(strfind(parameter_names{index_p}, 'centre')) ...
        |    ~isempty(strfind(parameter_names{index_p}, 'center')) ...
        |    ~isempty(strfind(parameter_names{index_p}, 'position'))
        if ~position1,     position1 = p(index_p);
        elseif ~position2, position2 = p(index_p);
        else               position2 = mean( [position2 p(index_p)] ); end % kind of center on all peaks
      end
    end
    if width1 == 0
      if ~position1, width1=1; else width1=abs(position1/2); end
    end
     if width2 == 0
      if ~position2, width2=1; else width2=abs(position2/2); end
    end
    signal.Axes   =  { linspace(position1-3*width1,position1+3*width1, 50) ...
                       linspace(position2-3*width2,position2+3*width2, 150) };
    signal.Values = evaluate(signal.Guess, signal.Axes{:});
  elseif nargin == 1 && ischar(p) && strcmp(p, 'plot')
    signal = feval(mfilename);
    if signal.Dimension == 1
      plot(signal.Axes{1}, signal.Values);
    elseif signal.Dimension == 2
      h=surf(signal.Axes{2}, signal.Axes{1}, signal.Values);set(h,'Edgecolor','none');
      xlabel('columns (2nd axis)'); ylabel('rows (1st axis)');
    end
    title(mfilename);
  elseif nargin == 0
    signal = feval(mfilename, [], linspace(-3,3, 150), linspace(-2,2, 100));
  else
    signal = identify;
  end

end
%                                        end of model main quad2d
% ------------------------------------------------------------------------------

% inline: evaluate: compute the model values
function signal = evaluate(p, x, y)
  if isempty(x) || isempty(y) || isempty(p), signal=[]; return; end
  % check if x,y are not mesh-like, then have to be computed on a mesh
  if isvector(x) && isvector(y)
    x=x(:); y=y(:)';
    [y,x]=meshgrid(y,x);
  end
  
  % HERE is the model evaluation <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  x0=p(2); y0=p(3); sx=p(4); sy=p(5);
  theta = p(6)*pi/180;  % deg -> rad
  a = cos(theta)^2/2/sx/sx + sin(theta)^2/2/sy/sy;
  b =-sin(2*theta)/4/sx/sx + sin(2*theta)/4/sy/sy;
  c = sin(theta)^2/2/sx/sx + cos(theta)^2/2/sy/sy;
  signal = (a*(x-x0).^2+2*b*(x-x0).*(y-y0)+c*(y-y0).^2);
  signal = p(1) *signal + p(7);

  signal = reshape(signal, size(x));
end

% inline: identify: return a structure which identifies the model
function signal =identify()
  % HERE are the parameter names
  parameter_names = {  'Amplitude' 'Centre_X' 'Center_Y' 'Curvature_X' 'Curvature_Y' 'Angle' 'Background' };
  %
  signal.Type           = 'iFit fitting function';
  signal.Name           = [ '2D Quadratic function with tilt angle (2D) [' mfilename ']' ];
  signal.Parameters     = parameter_names;
  signal.Dimension      = 2;         % dimensionality of input space (axes) and result
  signal.Guess          = [];        % default parameters
  signal.Axes           = {};        % the axes used to get the values
  signal.Values         = [];        % default model values=f(p)
  signal.function       = mfilename;
end

% inline: guess: guess some starting parameter values and return a structure
function info=guess(x,y,signal)
  info       = identify;  % create identification structure
  info.Axes  = { x y };
  % fill guessed information
  signal = signal(:); x=x(:); y=y(:);
  info.Guess = [ max(signal)-min(signal) mean(x) mean(y) std(x) std(y) 0 min(signal) ];
  info.Values= evaluate(info.Guess, info.Axes{:});
end

