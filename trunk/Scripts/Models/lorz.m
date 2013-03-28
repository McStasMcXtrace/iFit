function y=lorz(varargin)
% y = lorz(p, x, [y]) : Lorentzian
%
%   iFunc/lorz Lorentzian fitting function
%     y = p(1) ./ (1+ ((x-p(2)).^2/p(3)^2) ) + p(4);
%
% input:  p: Lorentzian model parameters (double)
%            p = [ Amplitude Centre HalfWidth BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=lorz([1 0 1 1], -10:10); or plot(lorz)
%
% Version: $Revision$
% See also iFunc, iFunc/fits, iFunc/plot

y.Name       = [ 'Lorentzian (1D) [' mfilename ']' ];
y.Parameters = {'Amplitude','Centre','HalfWidth','Background'};
y.Description= 'Single 1D Lorentzian model';
y.Expression = @(p,x) p(1) ./ (1+ ((x-p(2)).^2/p(3)^2) ) + p(4);
y.Guess     = @(x,signal) [ NaN ...
                            sum(signal(:).*x(:))/sum(signal(:)) ...
                            sqrt(abs(sum(x(:).*x(:).*signal(:))/sum(signal(:)) - sum(signal(:).*x(:))/sum(signal(:))*sum(signal(:).*x(:))/sum(signal(:)))) ...
                            NaN ];

y = iFunc(y);

if length(varargin)
  y = y(varargin{:});
end
