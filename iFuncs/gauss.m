function y=gauss(p, x, y)
% y = gauss(p, x, [y]) : Gaussian
%
%   iFunc/gauss Gaussian fitting function
%     y = p(1)*exp(-0.5*((x-p(2))/p(3)).^2) + p(4);
%   The function called with a char argument performs specific actions.
%   You may create new fit functions with the 'ifitmakefunc' tool.
%
% input:  p: Gaussian model parameters (double)
%            p = [ Amplitude Centre HalfWidth BackGround ]
%          or action e.g. 'identify', 'guess', 'plot' (char)
%         x: axis (double)
%         y: when values are given, a guess of the parameters is performed (double)
% output: y: model value or information structure (guess, identify)
% ex:     y=gauss([1 0 1 1], -10:10); or y=gauss('identify') or p=gauss('guess',x,y);
%
% Version: $Revision: 1.11 $
% See also iData, ifitmakefunc

% 1D function template:
% Please retain the function definition structure as defined below
% in most cases, just fill-in the information when HERE is indicated

  if nargin >= 2 && isnumeric(p) && ~isempty(p) && isnumeric(x) && ~isempty(x)
  %   evaluate: model(p,x, ...)
    y = evaluate(p, x);
  elseif nargin == 3 && isnumeric(x) && isnumeric(y) && ~isempty(x) && ~isempty(y)
  %   guess: model('guess', x,y)
  %   guess: model(p,       x,y)
    y = guess(x,y);
  elseif nargin == 2 && isnumeric(p) && isnumeric(x) && numel(p) == numel(x)
  %   guess: model(x,y) with numel(x)==numel(y)
    y = guess(p,x);
  elseif nargin == 2 && isnumeric(p) && ~isempty(p) && isempty(x)
  %   evaluate: model(p,[])
    y = feval(mfilename, p);
  elseif nargin == 2 && isempty(p) && isnumeric(x) && ~isempty(x)
  %   identify: model([],x)
    y = identify; x=x(:);
    % HERE default parameters when only axes are given <<<<<<<<<<<<<<<<<<<<<<<<<
    y.Guess  = [1 mean(x) std(x)/2 .1];
    y.Axes   = { x };
    y.Values = evaluate(y.Guess, y.Axes{:});
  elseif nargin == 1 && isnumeric(p) && ~isempty(p) 
  %   identify: model(p)
    y = identify;
    y.Guess  = p;
    % HERE default axes to represent the model when parameters are given <<<<<<<
    y.Axes   =  { linspace(p(2)-3*p(3),p(2)+3*p(3), 100) };
    y.Values = evaluate(y.Guess, y.Axes{:});
  elseif nargin == 1 && ischar(p) && strcmp(p, 'plot') % only works for 1D
    y = feval(mfilename, [], linspace(-2,2, 100));
    if y.Dimension == 1
      plot(y.Axes{1}, y.Values);
    elseif y.Dimension == 2
      surf(y.Axes{1}, y.Axes{2}, y.Values);
    end
    title(mfilename);
  elseif nargin == 0
    y = feval(mfilename, [], linspace(-2,2, 100));
  else
    y = identify;
  end

end
% end of model main
% ------------------------------------------------------------------------------

% inline: evaluate: compute the model values
function y = evaluate(p, x)
  sx = size(x); x=x(:);
  if isempty(x) | isempty(p), y=[]; return; end
  
  % HERE is the model evaluation <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  y = p(1)*exp(-0.5*((x-p(2))/p(3)).^2) + p(4);
  
  y = reshape(y, sx);
end

% inline: identify: return a structure which identifies the model
function y =identify()
  % HERE are the parameter names
  parameter_names = {'Amplitude','Centre','HalfWidth','Background'};
  %
  y.Type           = 'iFit fitting function';
  y.Name           = [ 'Gaussian (1D) [' mfilename ']' ];
  y.Parameters     = parameter_names;
  y.Dimension      = 1;         % dimensionality of input space (axes) and result
  y.Guess          = [];        % default parameters
  y.Axes           = {};        % the axes used to get the values
  y.Values         = [];        % default model values=f(p)
  y.function       = mfilename;
end

% inline: guess: guess some starting parameter values and return a structure
function info=guess(x,y)
  info       = identify;  % create identification structure
  info.Axes  = { x };
  % fill guessed information
  info.Guess = iFuncs_private_guess(x(:), y(:), info.Parameters);
  info.Values= evaluate(info.Guess, info.Axes{:});
end

