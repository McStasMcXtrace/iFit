function y=gauss(varargin)
% y = gauss(p, x, [y]) : Gaussian
%
%   iFunc/gauss Gaussian fitting function
%     y = p(1)*exp(-0.5*((x-p(2))/p(3)).^2) + p(4);
%
% input:  p: Gaussian model parameters (double)
%            p = [ Amplitude Centre HalfWidth BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=gauss([1 0 1 1], -10:10); or plot(gauss)
%
% Version: $Revision: 1.14 $
% See also iFunc, iFunc/fits, iFunc/plot

y.Name      = [ 'Gaussian (1D) [' mfilename ']' ];
y.Description='Single 1D Gaussian model';
y.Parameters={'Amplitude','Centre','HalfWidth','Background'};
y.Expression= @(p,x) p(1)*exp(-0.5*((x-p(2))/p(3)).^2) + p(4);
y.Guess     = @(x,signal) [ NaN ...
                            sum(signal(:).*x(:))/sum(signal(:)) ...
                            sqrt(abs(sum(x(:).*x(:).*signal(:))/sum(signal(:)) - sum(signal(:).*x(:))/sum(signal(:))*sum(signal(:).*x(:))/sum(signal(:)))) ...
                            NaN ];
y = iFunc(y);

if length(varargin)
  y = y(varargin{:});
end

