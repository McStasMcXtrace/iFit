function y=poisson(varargin)
% y = poisson(p, x, [y]) : Poisson distribution function. 
%
%   iFunc/poisson Poisson distribution function. 
%     y  = p(3)+ p(1)*(exp(-p(2)).*(p(2).^x))./factorial(floor(x));
%     valid with integer axis values
%
% poisson(centre)         creates a model with a specified centre
% poisson([ parameters ]) creates a model with specified model parameters
%
% input:  p: Poisson model parameters (double)
%            p = [ Amplitude Center BackGround ]
%          or action e.g. 'identify', 'guess', 'plot' (char)
%         x: axis (integers)
%         y: when values are given, a guess of the parameters is performed (double)
% output: y: model value or information structure (guess, identify)
% ex:     y=poisson([1 0 1 1], -10:10); or plot(poisson);
%
% Version: $Revision$
% See also iFunc, iFunc/fits, iFunc/plot

y.Name      = [ 'Poisson distribution function (1D) [' mfilename ']' ];
y.Parameters={'Amplitude','Center','Background'};
y.Description='Poisson distribution function. Only valid with integer axis values';
y.Expression= @(p,x) p(3)+ p(1)*(exp(-p(2)).*(p(2).^floor(abs(x))))./factorial(floor(abs(x)));
y.Guess     = @(x,y) [ (max(y(:))-min(y(:)))/2 mean(abs(x(:))) std(x(:)) ];
y.Dimension =1;

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 varargin{1} 0]};
  end
  y.ParameterValues = varargin{1};
elseif length(varargin) > 1
  y = y(varargin{:});
end

