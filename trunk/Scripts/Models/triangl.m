function y=triangl(varargin)
% y = triangl(p, x, [y]) : Triangular
%
%   iFunc/triangl Triangular fitting function
%     y=(p(3)-sign(x-p(2)).*(x-p(2)))/p(3)^2;
%     and y is set to the background outside the full width.
%
% triangl(width)          creates a model with a specified width
% triangl([ parameters ]) creates a model with specified model parameters
%
% input:  p: Triangular model parameters (double)
%            p = [ Amplitude Centre HalfWidth BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=triangl([1 0 1 1], -10:10); or plot(triangl);
%
% Version: $Revision$
% See also iFunc, iFunc/fits, iFunc/plot, heavisde, tophat

y.Name      = [ 'Triangular function (1D) [' mfilename ']' ];
y.Parameters={'Amplitude','Centre','HalfWidth','Background'};
y.Description='Triangular function';
% MZ <mzinkin@sghms.ac.uk>
y.Expression= @(p,x) ((p(3)-sign(x-p(2)).*(x-p(2)))/p(3)^2)*p(1).*(abs(x-p(2)) < p(3))+p(4);

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

