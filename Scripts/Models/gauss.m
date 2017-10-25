function y=gauss(varargin)
% y = gauss(p, x, [y]) : Gaussian
%
%   iFunc/gauss Gaussian fitting function
%     y = p(1)*exp(-0.5*((x-p(2))/p(3)).^2) + p(4);
%
%   The HalfWidth parameter is the Gaussian square root variance (Sigma). 
%   The 'true' half width is thus 1.177*HalfWidth.
%
% gauss(width)          creates a model with a specified width(sigma)
% gauss([ parameters ]) creates a model with specified model parameters
%
% Reference: http://en.wikipedia.org/wiki/Gaussian_function
%
% input:  p: Gaussian model parameters (double)
%            p = [ Amplitude Centre HalfWidth BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=gauss([1 0 1 1], -10:10); or plot(gauss)
%
% Version: $Date$
% See also iFunc, iFunc/fits, iFunc/plot, gauss1, lorz
% (c) E.Farhi, ILL. License: EUPL.

y.Name      = [ 'Gaussian (1D) [' mfilename ']' ];
y.Description='1D Gaussian model';
y.Parameters={'Amplitude','Centre','HalfWidth half Sigma RMS','Background'};
y.Expression= @(p,x) p(1)*exp(-0.5*((x-p(2))/p(3)).^2) + p(4);

% moments of distributions
m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));

y.Guess     = @(x,s) [ NaN m1(x, s-min(s(:))) m2(x, s-min(s(:))) NaN ];

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 0 varargin{1} 0]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

