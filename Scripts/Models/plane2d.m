function signal=plane2d(varargin)
% signal = plane2d(p, x, y, {signal}) : Planar function
%
%   iFunc/plane2d Planar function (fit 2D function/model)
%       signal = p(1)*x+p(2)*y+p(3)
%
% plane2d([s1 s2])        creates a model with specified slopes
% plane2d([ parameters ]) creates a model with specified model parameters
%
% input:  p: plane2d model parameters (double array)
%            p = [  'Slope_X' 'Slope_Y' 'Background' ]
%          or 'guess'
%         x: axis along rows    (double)
%         y: axis along columns (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value
% ex:     signal=plane2d([1 2 .5 .2 .3 30 .2], -2:.1:2, -3:.1:3); or plot(plane2d);
%
% Version: $Revision$
% See also iData, iFunc/fits, iFunc/plot, gauss

signal.Name           = [ 'Planar function (2D) [' mfilename ']' ];
signal.Parameters     = { 'Slope_X' 'Slope_Y' 'Background' };
signal.Description    = '2D Planar function';
signal.Dimension      = 2;         % dimensionality of input space (axes) and result
signal.Guess          = [];        % default parameters
signal.Expression     = @(p,x,y) p(1)*x+p(2)*y+p(3);

signal=iFunc(signal);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ varargin{1} varargin{1} 0]};
  elseif length(varargin{1}) == 2
    varargin = {[ varargin{1} 0]};
  end
  y.ParameterValues = varargin{1};
elseif length(varargin) > 1
  y = y(varargin{:});
end

