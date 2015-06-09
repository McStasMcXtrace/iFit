function y=sine(varargin)
% y = sine(p, x, [y]) : Sine function. 
%
%   iFunc/sine Sine function. 
%     y  = p(4) + p(1)*sin((x - p(2))/p(3))
%
% sine(period)         creates a model with specified period
% sine([ parameters ]) creates a model with specified model parameters
%
% input:  p: Sine model parameters (double)
%            p = [ Amplitude Phase_Shift Period BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=sine([1 0 1 1], -10:10); or plot(sine);
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot

y.Name           = [ 'Sine function (1D) [' mfilename ']' ];
y.Parameters     = {'Amplitude','Phase_Shift','Period','Background'};
y.Dimension      = 1;
y.Description    = 'Sine function';
y.Expression     = @(p,x) p(4) + p(1)*sin((x - p(2))/p(3));
y.Guess          =  @(x,y) [ max(y(:))-min(y(:)) mean(x(:)) std(x(:))/4 min(y(:)) ];
y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 0 varargin{1} 0]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

