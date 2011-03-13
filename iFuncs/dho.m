function y=dho(p, x, y)
% y = dho(p, x, [y]) : Damped harmonic oscillator
%
%   iFunc/dho Damped harmonic oscillator fitting function, including Bose factor
%   the function called with a char argument performs specific actions.
%
% input:  p: Damped harmonic oscillator model parameters (double)
%            p = [ Amplitude Centre HalfWidth BackGround Temperature(in x units)]
%          or action e.g. 'identify', 'guess' (char)
%         x: axis (double)
%         y: when values are given, a guess of the parameters is performed (double)
% output: y: model value or information structure (guess, identify)
% ex:     y=dho([1 0 1 1], -10:10); or y=dho('identify') or p=dho('guess',x,y);

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
    y.Guess  = [1 mean(x) std(x)/2 .1 1];
    y.Axes   = { x };
    y.Values = evaluate(y.Guess, y.Axes{:});
  elseif nargin == 1 && isnumeric(p) && ~isempty(p) 
  %   identify: model(p)
    y = identify;
    y.Guess  = p;
    % HERE default axes to represent the model when parameters are given <<<<<<<
    y.Axes   =  { linspace(p(2)-3*p(3),p(2)+3*p(3), 100) };
    y.Values = evaluate(y.Guess, y.Axes{:});
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
  % DFM 19.10.95 rev EF 27.06.97
  omega=x(:);
  amp=p(1);
  omega_0=p(2); 
  gamma=p(3);	
  T=p(5);

  n_omega=1./(exp(abs(omega)/T)-1);	% bose

  s_omega = 0*x;
  as_omega = 0*x;

  i = find(omega > 0);

  s_omega(i)=gamma*omega(i)*omega_0^2.*(1+n_omega(i))./((omega(i).^2-omega_0^2).^2+(gamma*omega(i)).^2); % dho stokes

  omega = -omega;

  i = find(omega > 0);

  as_omega(i) =gamma*omega(i)*omega_0^2.*(n_omega(i))./((omega(i).^2-omega_0^2).^2 + (gamma*omega(i)).^2); % dho antistokes

  y=amp*(s_omega+as_omega)+p(4);
  
  y = reshape(y, sx);
end

% inline: identify: return a structure which identifies the model
function y =identify()
  % HERE are the parameter names
  parameter_names = {'Amplitude','Centre','HalfWidth','Background','T ("x" unit)'};
  %
  y.Type           = 'iFit fitting function';
  y.Name           = [ 'Damped harmonic oscillator (1D) [' mfilename ']' ];
  y.Parameters     = parameter_names;
  y.Dimension      = 1;         % dimensionality of input space (axes) and result
  y.Guess          = [];        % default parameters
  y.Axes           = {};        % the axes used to get the values
  y.Values         = [];        % default model values=f(p)
end

% inline: guess: guess some starting parameter values and return a structure
function info=guess(x,y)
  info       = identify;  % create identification structure
  info.Axes  = { x };
  % fill guessed information
  info.Guess = iFuncs_private_guess(x(:), y(:), info.Parameters);
  info.Guess(end) = 1;
  info.Values= evaluate(info.Guess, info.Axes{:});
end

