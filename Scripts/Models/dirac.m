function y=dirac(varargin)
% y = dirac(p, x, [y]) : Dirac
%
%   iFunc/dirac Dirac fitting function
%      y(x==p(2)) = p(1)
%
% dirac(centre)         creates a model with a line at centre
% dirac([ parameters ]) creates a model with specified model parameters
%
% input:  p: Dirac model parameters (double)
%            p = [ Amplitude Centre ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=dirac([1 0], -10:10); or plot(dirac)
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot
% (c) E.Farhi, ILL. License: EUPL.

y.Name       = [ 'Dirac (1D) [' mfilename ']' ];
y.Description='Dirac peak fitting function';
y.Parameters = {'Amplitude','Centre'};
y.Expression = @(p,x) p(1)*(abs(x - p(2)) == min(abs(x(:) - p(2))));
y.Dimension  = 1;
% moments of distributions
m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));

y.Guess      = @(x,signal) [ NaN m1(x, signal-min(signal(:))) ];

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 varargin{1} ]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end


