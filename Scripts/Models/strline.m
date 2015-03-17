function y=strline(varargin)
% y = strline(p, x, [y]) : Straight line
%
%   iFunc/strline Straight line fitting function
%     y=p(2)+p(1)*x;
%     The p(1)=Gradient parameter is the slope of the straight line.
%
% strline(slope)          creates a model with specified slope/gradient
% strline([ parameters ]) creates a model with specified model parameters
%
% input:  p: Straight line model parameters (double)
%            p = [ Gradient BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=strline([1 1], -10:10); or plot(strline);
%
% Version: $Revision$
% See also iFunc, iFunc/fits, iFunc/plot, quadline, plane2d

y.Guess     = @(x,y) [ polyfit(x(:), y(:), 1) ];
y.Name      = [ 'Straight-line (1D) [' mfilename ']' ];
y.Parameters={'Gradient slope','Constant'};
y.Description='Straight line';
y.Expression= @(p,x) p(2)+p(1)*x;

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ varargin{1} 0 ]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

