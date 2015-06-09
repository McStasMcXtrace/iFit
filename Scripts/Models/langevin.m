function y=langevin(varargin)
% y = langevin(p, x, [y]) : Langevin function for magnetic polarization
%
%   iFunc/langevin Langevin function for magnetic polarization
%     x = (x-p(2))/p(3); y  = p(4)+ p(1) .* ( coth(x) - 1./x );
%
% langevin(centre)         creates a model with a specified centre/threshold
% langevin([ parameters ]) creates a model with specified model parameters
%
% input:  p: Langevin model parameters (double)
%            p = [ Amplitude Center Width BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=langevin([1 0 1 1], -10:10); or plot(langevin);
%
% Version: $Date$
% See also iFunc, iFunc/fits, iFunc/plot

y.Name      = [ 'Langevin function for magnetic polarization (1D) [' mfilename ']' ];
y.Parameters={'Amplitude','Centre','Width','Background'};
y.Description='Langevin function for magnetic polarization';
y.Expression= @(p,x) p(4)+ p(1) .* ( coth((x-p(2))/p(3)) - p(3)./(x-p(2)) );
y.Guess     = @(x,y) [ (max(y(:))-min(y(:)))/2 mean(x(:)) std(x(:)) min(y(:)) ];

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 varargin{1} varargin{1}/4 0]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

