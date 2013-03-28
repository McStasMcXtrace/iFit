function y=bilorz(varargin)
% y = bilorz(p, x, [y]) : Bi-Lorentzian
%
%   iFunc/bilorz Bi-Lorentzian fitting function
%     y = p(1)*exp(-0.5*((x-p(2))/s).^2) + p(5);
%   where s = p(3) for x < p(2) and s = p(4) for x > p(2).
%
% input:  p: Bi-Lorentzian model parameters (double)
%            p = [ Amplitude Centre HalfWidth1 HalfWidth2 BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=bilorz([1 0 1 1], -10:10); or plot(bilorz)
%
% Version: $Revision$
% See also iFunc, iFunc/fits, iFunc/plot

y.Name      = [ 'Bi-Lorentzian (1D) [' mfilename ']' ];
y.Description='Bi-Lorentzian/asymmetric fitting function';
y.Parameters= { 'Amplitude', 'Centre', 'HalfWidth1', 'HalfWidth2', 'BackGround' };
y.Expression = @(p,x) p(1)*exp(-0.5*((x-p(2))./ (p(3)*(x < p(2)) + p(4) * (x >= p(2))) ).^2) + p(5);
y.Guess     = @(x,signal) [ NaN ...
                            sum(signal(:).*x(:))/sum(signal(:)) ...
                            sqrt(abs(sum(x(:).*x(:).*signal(:))/sum(signal(:)) - sum(signal(:).*x(:))/sum(signal(:))*sum(signal(:).*x(:))/sum(signal(:))))/1.5 ...
                            sqrt(abs(sum(x(:).*x(:).*signal(:))/sum(signal(:)) - sum(signal(:).*x(:))/sum(signal(:))*sum(signal(:).*x(:))/sum(signal(:))))*1.5 ...
                            NaN ];

y=iFunc(y);

if length(varargin)
  y = y(varargin{:});
end
