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
% ex:     y=poisson([1 1 0], -10:10); or plot(poisson);
%
% Version: $Date$
% See also iFunc, iFunc/fits, iFunc/plot
% (c) E.Farhi, ILL. License: EUPL.

y.Name      = [ 'Poisson distribution function (1D) [' mfilename ']' ];
y.Parameters={'Amplitude','Center','Background'};
y.Description='Poisson distribution function. Only valid with integer axis values';
y.Expression= @(p,x) p(3)+ p(1)*(exp(-p(2)).*(p(2).^floor(abs(x))))./factorial(floor(abs(x)));

% use ifthenelse anonymous function
% <https://blogs.mathworks.com/loren/2013/01/10/introduction-to-functional-programming-with-anonymous-functions-part-1/>
% iif( cond1, exec1, cond2, exec2, ...)
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
y.Guess     = @(x,y) iif(...
  ~isempty(y), @() [ (max(y(:))-min(y(:)))/2 mean(abs(x(:))) min(x(:)) ], ...
  true            , @() [1 1 0]);
  
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

