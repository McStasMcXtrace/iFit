function y=gauss1(varargin)
% y = gauss1(p, x, [y]) : Gaussian without background
%
%   iFunc/gauss1 Gaussian fitting function without background
%     sigma = p(3)/(2*sqrt(2*log(2)))=p(3)/2.3548
%     y = p(1)/(sigma*sqrt(2*pi))*exp(-0.5*((x-p(2))/sigma).^2);
%
%     The Amplitude is:       p(1)/(p(3)/(2*sqrt(2*log(2)))*sqrt(2*pi))
%
% gauss1(width)          creates a model with a specified width
% gauss1([ parameters ]) creates a model with specified model parameters
%
% Reference: http://en.wikipedia.org/wiki/Gaussian_function
%
% input:  p: Gaussian model parameters (double)
%            p = [ Intensity Centre HalfWidth ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=gauss1([1 0 1], -10:10); or plot(gauss1)
%
% Version: $Date$
% See also iFunc, iFunc/fits, iFunc/plot, lorz, gauss, lorz1
% (c) E.Farhi, ILL. License: EUPL.

y.Name      = [ 'Gaussian1 (1D) [' mfilename ']' ];
y.Description='1D Gaussian1 model';
y.Parameters={'Intensity','Centre','HalfWidth 2.3548*sigma'};
y.Expression= @(p,x) p(1)/(p(3)/(2*sqrt(2*log(2)))*sqrt(2*pi))*exp(-0.5*((x-p(2))/p(3)*(2*sqrt(2*log(2)))).^2);
y.Dimension  = 1;

% moments of distributions
m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));

% use ifthenelse anonymous function
% <https://blogs.mathworks.com/loren/2013/01/10/introduction-to-functional-programming-with-anonymous-functions-part-1/>
% iif( cond1, exec1, cond2, exec2, ...)
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
y.Guess     = @(x,s) iif(...
  ~isempty(s)&&numel(s)==numel(x), @()  [ NaN m1(x, s-min(s(:))) m2(x, s-min(s(:))) ], ...
  true            , @() [1 0 1]);

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 0 varargin{1} ]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

