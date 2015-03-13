function y=gauss(varargin)
% y = gauss(p, x, [y]) : Gaussian
%
%   iFunc/gauss Gaussian fitting function
%     y = p(1)*exp(-0.5*((x-p(2))/p(3)).^2) + p(4);
%
% gauss(width)          creates a model with a specified width
% gauss([ parameters ]) creates a model with specified model parameters
%
% input:  p: Gaussian model parameters (double)
%            p = [ Amplitude Centre HalfWidth BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=gauss([1 0 1 1], -10:10); or plot(gauss)
%
% Version: $Revision$
% See also iFunc, iFunc/fits, iFunc/plot

y.Name      = [ 'Gaussian (1D) [' mfilename ']' ];
y.Description='1D Gaussian model';
y.Parameters={'Amplitude','Centre','HalfWidth','Background'};
y.Expression= @(p,x) p(1)*exp(-0.5*((x-p(2))/p(3)).^2) + p(4);

% moments of distributions
m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));

y.Guess     = @(x,s) [ NaN m1(x, s-min(s(:))) m2(x, s-min(s(:))) NaN ];

y = iFunc(y);

%if nargin == 1 && isnumeric(varargin{1})
  %if length(varargin{1}) == 1
    %varargin = {[ 1 0 varargin{1} 0]};
  %end
  %y.ParameterValues = varargin{1};
%end

if length(varargin)
  y = y(varargin{:});
end

