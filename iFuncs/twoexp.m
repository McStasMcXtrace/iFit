function signal=twoexp(p, x, signal)
% signal = twoexp(p, x, {signal}) : two exponential decay functions
%
%   iFunc/twoexp two exponential decay functions (fit 1D function/model)
%     y=p(1)*exp(-x/p(2))+p(3)*exp(-x/p(4)) + p(5);
%   the function called with a char argument performs specific actions.
%
% input:  p: twoexp model parameters (double)
%            p = [   'Amplitude1' 'Tau1' 'Amplitude2' 'Tau2' 'Background' ] as a numerical array
%          or action e.g. 'identify', 'guess', 'plot' (char)
%         x: axis (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value or information structure ('guess', 'identify','plot')
% ex:     signal=twoexp([1 0 1 1], -10:10); or signal=twoexp('identify') or p=twoexp('guess',x,signal);
%
% Version: $Revision: 1.1 $
% See also iData, ifitmakefunc

% 1D function created by ifitmakefunc from template (iFit/iFuncs)
%   Please retain the function definition structure as defined below
%   in most cases, just fill-in the information when HERE is indicated

  if nargin >= 2 && isnumeric(p) && ~isempty(p) && isnumeric(x) && ~isempty(x)
  %   evaluate: model(p,x, ...)
    signal = evaluate(p, x);
  elseif nargin == 3 && isnumeric(x) && isnumeric(signal) && ~isempty(x) && ~isempty(signal)
  %   guess: model('guess', x,signal)
  %   guess: model(p,       x,signal)
    signal = guess(x,signal);
  elseif nargin == 2 && isnumeric(p) && isnumeric(x) && numel(p) == numel(x)
  %   guess: model(x,signal) with numel(x)==numel(signal)
    signal = guess(p,x);
  elseif nargin == 2 && isnumeric(p) && ~isempty(p) && isempty(x)
  %   evaluate: model(p,[])
    signal = feval(mfilename, p);
  elseif nargin == 2 && isempty(p) && isnumeric(x) && ~isempty(x)
  %   identify: model([],x)
    info = identify; x=x(:); signal=[];
    info.Guess  = [ 1 std(x)/2 0.5 std(x)*2 0.01 ];
    info.Axes   = { x };
    info.Values = evaluate(info.Guess, info.Axes{:});
    signal = info;
  elseif nargin == 1 && isnumeric(p) && ~isempty(p) 
  %   identify: model(p)
    signal = identify;
    signal.Guess  = p;
    width    = 0;
    position = 0;
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
    signal.Axes   =  { linspace(position-3*width,position+3*width, 100) };
    signal.Values = evaluate(signal.Guess, signal.Axes{:});
  elseif nargin == 1 && ischar(p) && strcmp(p, 'plot') % only works for 1D
    signal = feval(mfilename, [], linspace(0,2, 100));
    if signal.Dimension == 1
      plot(signal.Axes{1}, signal.Values);
    elseif signal.Dimension == 2
      surf(signal.Axes{1}, signal.Axes{2}, signal.Values);
    end
    title(mfilename);
  elseif nargin == 0
    signal = feval(mfilename, [], linspace(0,2, 100));
  else
    signal = identify;
  end

end
%                                        end of model main twoexp
% ------------------------------------------------------------------------------

% inline: evaluate: compute the model values
function signal = evaluate(p, x)
  sx = size(x); x=x(:);
  if isempty(x) | isempty(p), signal=[]; return; end
  
  % HERE is the model evaluation <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    signal = p(1)*exp(-x/p(2))+p(3)*exp(-x/p(4)) + p(5);
  
  signal = reshape(signal, sx);
end

% inline: identify: return a structure which identifies the model
function signal =identify()
  % HERE are the parameter names
  parameter_names = {  'Amplitude1' 'Tau1' 'Amplitude2' 'Tau2' 'Background'};
  %
  signal.Type           = 'iFit fitting function';
  signal.Name           = [ 'two exponential decay functions (1D) [' mfilename ']' ];
  signal.Parameters     = parameter_names;
  signal.Dimension      = 1;         % dimensionality of input space (axes) and result
  signal.Guess          = [ 1 2 0.5 3 0.01 ];        % default parameters
  signal.Axes           = {};        % the axes used to get the values
  signal.Values         = [];        % default model values=f(p)
  signal.function       = mfilename;
end

% inline: guess: guess some starting parameter values and return a structure
function info=guess(x,y)
  info       = identify;  % create identification structure
  info.Axes  = { x };
  % fill guessed information
  mny = min(y);
  if min(y) <= 0, bkg = min(y); else bkg = 0; end
  pe = polyfit(x,log(y - bkg - abs(bkg)*0.01),1); % p(1) is the highest degree parameter
  p(1) = exp(pe(2))/2; p(3) = p(1);
  p(2) = -1/pe(1); p(4) = p(2);
  p(5) = mny;
  info.Guess = p;
  info.Values= evaluate(info.Guess, info.Axes{:});
end

