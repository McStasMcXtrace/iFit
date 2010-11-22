function y=gauss(p, varargin)
% y = gauss(p, x) : Gaussian
%
%   iFunc/gauss Gaussian fitting function
%   the function called with a char argument performs specific actions.
%
% input:  p: Gaussian model parameters
%            p = [ Amplitude Centre HalfWidth BackGround ]
%         x: axis (double) or action e.g. 'identify' (char)
% output: y: model value
% ex:     y=gauss([1 0 1 1], -10:10); or y=gauss('identify') or p=gauss('guess',x,y);

% Please retain the function definition structure which should provide:
% * model evaluation
% * model identification
% * guessed parameters from input data

if nargin==2 && isnumeric(p)
  % use of model(x,y) triggers parameter guess
  x = varargin{1};
  if length(p) == length(x)
    y = feval(mfilename, 'guess', p, x);
    return
  end
  
  % use of model([], x) triggers default parameters to be used
  if isempty(p) 
    id=feval(mfilename, 'identify');
    p=id.Guess;
  end
  
  % evaluate model values
  y=p(1)*exp(-0.5*((x-p(2))/p(3)).^2) + p(4);
  
else
  % check for actions, default is identify
  action='';
  if nargin == 0, p=[]; end
  if ischar(p), action=p; end
  if isempty(action), action='identify'; end
  
  % the names of the model parameters
  parameter_names = {'Amplitude','Centre','HalfWidth','Background'};
  
  % some reasonable X axis to plot the function when nothing is given
  if nargin < 2
    % if parameters have been given use them to define a sensible X axis
    if ~isempty(p)
      x = linspace(p(2)-3*p(3),p(2)+3*p(3), 100);
    else
    % else use default
      x = linspace(-2,2,100);
    end
  end
  
  % when we provide the data, guess some parameters y =model(p,x)
  if nargin >= 3
    x = varargin{1};
    y = varargin{end};
    % use the Y values to guess...

    % make a peak search in 1D signal y=f(x) and assign parameters according to their names
    p = iFuncs_private_guess(x, y, parameter_names);
    y = [];
  else
    % default parameters when nothing is given as argument
    p = [1 mean(x) std(x)/2 .1];
  end
  
  % compute the model values with default, given or guessed parameters
  model = feval(mfilename, p, x);
  
  % some actions for this function to define return value
  switch action
  case 'identify'
    y.Type           = 'iFit fitting function';
    y.Name           = 'Gaussian (1D)';
    y.Parameters     = parameter_names;
    y.Dimension      = 1;         % dimensionality of input space (axes) and result
    y.Guess          = p;         % default or guessed parameters
    y.Values         = model;     % default or guessed model values=f(p)
    y.Axis           = { x };     % the axes used to get the values
  case 'guess'
    y = p;
  otherwise
    disp([ mfilename ': unknown action ' action ]);
    y = [];
  end
  % and we return 'y' 
end

