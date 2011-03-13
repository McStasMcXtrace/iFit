function y=voigt(p, x, y)
% y = voigt(p, x, [y]) : Voigt function
%
%   iFunc/voigt Voigt fitting function, including Bose factor
%   the function called with a char argument performs specific actions.
%
% input:  p: Voigt model parameters (double)
%            p = [ Amplitude Centre Width_Gauss Width_Lorz BackGround ]
%          or action e.g. 'identify', 'guess' (char)
%         x: axis (double)
%         y: when values are given, a guess of the parameters is performed (double)
% output: y: model value or information structure (guess, identify)
% ex:     y=voigt([1 0 1 1], -10:10); or y=voigt('identify') or p=voigt('guess',x,y);

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
    y.Guess  = [1 mean(x) std(x)/2 std(x)/2 .1 1];
    y.Axes   = { x };
    y.Values = evaluate(y.Guess, y.Axes{:});
  elseif nargin == 1 && isnumeric(p) && ~isempty(p) 
  %   identify: model(p)
    y = identify;
    y.Guess  = p;
    % HERE default axes to represent the model when parameters are given <<<<<<<
    y.Axes   =  { linspace(p(2)-3*(p(3)+p(4)),p(2)+3*(p(3)+p(4)), 100) };
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
  % MZ <mzinkin@sghms.ac.uk> adapted from DFM
  N = 16;
	b = -sqrt(log(2))/p(3);
	a = b*p(4);
	b = b*2*i;
	z = a + b*(x-p(2));

	M=2*N; M2=2*M; k=[-M+1:1:M-1]';
	L=sqrt(N/sqrt(2));
	tt=(L*tan(k*pi/M2)).^2;
	f=[0; exp(-tt).*(L^2+tt)];
	a=real(fft(fftshift(f)))/M2;
	a=flipud(a(2:N+1));
	l=L-z;
	Z=(L+z)./l;
	pp=polyval(a,Z);
	y=p(5)+p(1)*real(2*pp ./l.^2+(1/sqrt(pi))*ones(size(z)) ./l);
end

% inline: identify: return a structure which identifies the model
function y =identify()
  % HERE are the parameter names
  parameter_names = {'Amplitude','Centre','Width_Gauss',' Width_Lorz', 'Background'};
  %
  y.Type           = 'iFit fitting function';
  y.Name           = [ 'Voigt (1D) [' mfilename ']' ];
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
  info.Guess(4) = info.Guess(3);
  info.Values= evaluate(info.Guess, info.Axes{:});
end

