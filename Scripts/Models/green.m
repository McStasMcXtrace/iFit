function y=green(varargin)
% y = green(p, x, [y]) : Green function
%
%   iFunc/green Green fitting function
%     y = (p(1)*p(3)*p(2)^2 ) ./ ( (p(2)^2 - x.^2).^2 + (x*p(3)).^2) + p(4);
%
% green(energy)         creates a model with a line at specified energy (centre)
% green([ parameters ]) creates a model with specified model parameters
%
% input:  p: Green model parameters (double)
%            p = [ Amplitude Center HalfWidth BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=green([1 0 1 0], -10:10); or plot(green);
%
% Version: $Date$
% See also iFunc, iFunc/fits, iFunc/plot
% $

y.Name      = [ 'Green function (1D) [' mfilename ']' ];
y.Parameters={'Amplitude','Centre','HalfWidth','Background'};
y.Description='Green function model';
y.Expression= @(p,x) (p(1)*p(3)*p(2)^2 ) ./ ( (p(2)^2 - x.^2).^2 + (x*p(3)).^2) + p(4);
y.Dimension  = 1;

% moments of distributions
m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));

% use ifthenelse anonymous function
% <https://blogs.mathworks.com/loren/2013/01/10/introduction-to-functional-programming-with-anonymous-functions-part-1/>
% iif( cond1, exec1, cond2, exec2, ...)
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
y.Guess     = @(x,s) iif(...
  ~isempty(s)&&numel(s)==numel(x), @()  [ NaN m1(x, s-min(s(:))) m2(x, s-min(s(:))) NaN ], ...
  true            , @() [1 0 1 0]);

y = iFunc(y);


if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 varargin{1} varargin{1}/4 0]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end


