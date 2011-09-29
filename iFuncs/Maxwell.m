function y=Maxwell(p, x, y)
% y = Maxwell(p, x, [y]) : Maxwellian
%
%   iFunc/Maxwell Maxwellian fitting function
%
%   The function called with a char argument performs specific actions.
%   You may create new fit functions with the 'ifitmakefunc' tool.
%
% input:  p: Maxwellian model parameters (double)
%            p = [ Amplitude Centre HalfWidth BackGround ]
%          or action e.g. 'identify', 'guess', 'plot' (char)
%         x: axis (double)
%         y: when values are given, a guess of the parameters is performed (double)
% output: y: model value or information structure (guess, identify)
% ex:     y=Maxwell([1 0 1 1], -10:10); or y=Maxwell('identify') or p=Maxwell('guess',x,y);
%
% Version: $Revision: 1.1 $
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
    y.Guess  = [213 1.46e13 83 3.26e12 26 1.2e13 ];
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
  HBAR    =1.05459E-34;
  MNEUTRON=1.67492E-27;
  k       = 1.38066e-23;
  lambda=x;
  T1=p(1);
  I1=p(2);
  T2=p(3);
  I2=p(4);
  T3=p(5);
  I3=p(6);
    if (T1>0)
      lambda0  = 1.0e10*sqrt(HBAR*HBAR*4.0*pi*pi/2.0/MNEUTRON/k/T1);
      lambda02 = lambda0*lambda0;	   
      L2P      = 2*lambda02*lambda02;
    else lambda0=0;
    end
    if (T2>0)
      lambda0b  = 1.0e10*sqrt(HBAR*HBAR*4.0*pi*pi/2.0/MNEUTRON/k/T2);
      lambda02b = lambda0b*lambda0b;	   
      L2Pb      = 2*lambda02b*lambda02b;
    else lambda0b=0;
    end
    if (T3>0)
      lambda0c  = 1.0e10*sqrt(HBAR*HBAR*4.0*pi*pi/2.0/MNEUTRON/k/T3);
      lambda02c = lambda0c*lambda0c;	   
      L2Pc      = 2*lambda02c*lambda02c;
    else lambda0c=0;
    end
    lambda2=lambda .*lambda;
    lambda5=lambda2.*lambda2.*lambda;
    Maxwell=I1*zeros(size(x));
    if (T1 > 0)
      Maxwell= I1*L2P./lambda5.*exp(-lambda02./lambda2);
      if ((T2 > 0) & (I1 ~= 0))
        if (I2 == 0), I2 = I1; end
        Maxwell = Maxwell+ (I2).*L2Pb./lambda5.*exp(-lambda02b./lambda2);
      end
      if ((T3 > 0) & (I1 ~= 0))
        if (I3 == 0), I3 = I1; end
        Maxwell = Maxwell+ (I3).*L2Pc./lambda5.*exp(-lambda02c./lambda2);
       end
    end
    y=Maxwell;
  
  y = reshape(y, sx);
end

% inline: identify: return a structure which identifies the model
function y =identify()
  % HERE are the parameter names
  parameter_names = {'T1','I1','T2','I2','T3','I3'};
  %
  y.Type           = 'iFit fitting function';
  y.Name           = [ 'Maxwellian (1D) [' mfilename ']' ];
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
  % compute first and second moment
  x = x(:); y=y(:);
  
  info.Values= evaluate(info.Guess, info.Axes{:});
end

