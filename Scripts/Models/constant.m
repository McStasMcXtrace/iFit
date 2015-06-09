function y=constant(varargin)
% y = constant(p) : Constant
%
%   iFunc/constant a constant/background/single value
%     y=p(1);
%
% constant(cte)          creates a model with specified constant
%
% input:  p: Constant parameter (double)
%            p = [ Value ]
%          or 'guess' 
%          or name of a new parameter (string)
% output: y: model value
% ex:     y=constant('Temperature');
%
% Version: $Date$
% See also iFunc, iFunc/fits, iFunc/plot, quadline, plane2d

if length(varargin)==1 && ischar(varargin{1}) ...
        && ~any(strcmp(varargin{1}, {'guess','identify'})) && ~isempty(varargin{1})
  y.Name      = varargin{1};
  y.Parameters= { y.Name };
  varargin(1) = [];
else
  y.Name      = [ 'Constant (0D) [' mfilename ']' ];
  y.Parameters= {'Constant'};
end

y.Guess     = [];
y.Description='Constant/User value';
y.Expression= @(p,x) p(1);

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  y.ParameterValues = varargin{1};
else
  y = y(varargin{:});
end

