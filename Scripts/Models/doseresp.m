function y=doseresp(varargin)
% y = doseresp(p, x, [y]) : Dose-response curve with variable Hill slope
%
%   iFunc/doseresp Dose-response curve with variable Hill slope (sigmoid)
%     y  = p(4)+ p(1) ./ (1+10.^((p(2)-x).*p(3)));
%   This is a sigmoid S-shaped curve, aka logistic.
%
% doseresp(threshold)      creates a model with a specified threshold
% doseresp([ parameters ]) creates a model with specified model parameters
%
% input:  p: Dose Response model parameters (double)
%            p = [ Amplitude Center Slope BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=doseresp([1 0 1 1], -10:10); or plot(doseresp)
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, sigmoid

y.Name           = [ 'Dose-response (sigmoid) (1D) [' mfilename ']' ];
y.Description    = 'sigmoid S-shaped curve, aka logistic, aka dose response';
y.Parameters     = {'Amplitude','Center threshold','Slope','Background'};
y.Expression = @(p,x) p(4)+ p(1) ./ (1+10.^((p(2)-x).*p(3)));
y.Dimension      = 1;
y.Guess          =  @(x,y) [ max(y(:))-min(y(:)) mean(x(:)) (max(y(:))-min(y(:)))/std(x(:)) min(y(:)) ];
y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 varargin{1} 1 0 ]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end


