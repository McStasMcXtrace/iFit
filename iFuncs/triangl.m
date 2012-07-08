function y=triangl(varargin)
% y = triangl(p, x, [y]) : Triangular
%
%   iFunc/triangl Triangular fitting function
%     y=(p(3)-sign(x-p(2)).*(x-p(2)))/p(3)^2;
%     and y is set to the background outside the full width.
%
% input:  p: Triangular model parameters (double)
%            p = [ Amplitude Centre HalfWidth BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=triangl([1 0 1 1], -10:10); or plot(triangl);
%
% Version: $Revision: 1.6 $
% See also iFunc, iFunc/fits, iFunc/plot, heavisde, tophat

y.Name      = [ 'Triangular function (1D) [' mfilename ']' ];
y.Parameters={'Amplitude','Centre','HalfWidth','Background'};
y.Description='Triangular function';
% MZ <mzinkin@sghms.ac.uk>
y.Expression= @(p,x) ((p(3)-sign(x-p(2)).*(x-p(2)))/p(3)^2)*p(1).*(abs(x-p(2)) < p(3))+p(4);
y.Guess     = @(x,signal) [ NaN ...
                            sum(signal(:).*x(:))/sum(signal(:)) ...
                            sqrt(abs(sum(x(:).*x(:).*signal(:))/sum(signal(:)) - sum(signal(:).*x(:))/sum(signal(:))*sum(signal(:).*x(:))/sum(signal(:)))) ...
                            NaN ];
                            
y = iFunc(y);

if length(varargin)
  y = y(varargin{:});
end

