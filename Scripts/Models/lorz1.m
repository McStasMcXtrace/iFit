function y=lorz1(varargin)
% y = lorz1(p, x, [y]) : Lorentzian without background
%
%   iFunc/lorz1 Lorentzian fitting function without background
%     y = p(1)*p(3)/2/pi ./ ((x-p(2)).^2 +(p(3)/2).^2);
%
%   The Amplitude is 2/pi/p(3)
%
% lorz1(width)          creates a model with a specified width
% lorz1([ parameters ]) creates a model with specified model parameters
%
% Reference: http://en.wikipedia.org/wiki/Lorentzian_function
%
% input:  p: Lorentzian model parameters (double)
%            p = [ Intensity Centre HalfWidth ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=lorz1([1 0 1], -10:10); or plot(lorz1)
%
% Version: $Date$
% See also iFunc, iFunc/fits, iFunc/plot, gauss, lorz
% (c) E.Farhi, ILL. License: EUPL.

y.Name       = [ 'Lorentzian1 (1D) [' mfilename ']' ];
y.Parameters = {'Intensity','Centre','HalfWidth'};
y.Description= '1D Lorentzian1 model';
y.Expression = @(p,x) p(1)*p(3)/2/pi ./ ((x-p(2)).^2 +(p(3)/2).^2);

% moments of distributions
m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));

y.Guess     = @(x,s) [ NaN m1(x, s-min(s(:))) m2(x, s-min(s(:))) ];

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 0 varargin{1}]};
  end
  y.ParameterValues = varargin{1};
elseif length(varargin) > 1
  y = y(varargin{:});
end
