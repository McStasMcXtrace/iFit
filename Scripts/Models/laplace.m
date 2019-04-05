function y=laplace(varargin)
% y = laplace(p, x, [y]) : Laplace distribution function. 
%
%   iFunc/laplace Laplace distribution function. 
%     y  = p(4)+ p(1)/2/p(3) .* exp( abs(x - p(2))/p(3) );
%
% laplace(width)          creates a model with a specified width
% laplace([ parameters ]) creates a model with specified model parameters
%
% input:  p: Laplace model parameters (double)
%            p = [ Amplitude Center Width BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=laplace([1 0 1 0], -10:10); or plot(laplace);
%
% Version: $Date$
% See also iFunc, iFunc/fits, iFunc/plot
% $

y.Name      = [ 'Laplace distribution function (1D) [' mfilename ']' ];
y.Parameters={'Amplitude','Centre','Width','Background'};
y.Description='Laplace distribution function';
y.Expression= @(p,x) p(4)+ p(1)/2/p(3) .* exp( - abs(x - p(2))/p(3) );
y.Dimension  = 1;

% use ifthenelse anonymous function
% <https://blogs.mathworks.com/loren/2013/01/10/introduction-to-functional-programming-with-anonymous-functions-part-1/>
% iif( cond1, exec1, cond2, exec2, ...)
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
y.Guess     = @(x,y) iif(...
  ~isempty(y), @()  [ (max(y(:))-min(y(:)))/2 mean(x(:)) std(x(:)) min(y(:)) ], ...
  true            , @() [1 0 1 0]);

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 0 varargin{1} 0]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

