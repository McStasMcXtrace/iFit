function y=bose(varargin)
% y = bose(p, x, [y]) : Bose
%
%   iFunc/bose Bose fitting function
%     n(x) = 1 ./ (exp(x/p(1)) - 1);
%     y = n(x)+1 for x > 0
%     y = n(x)   for x < 0
%         
% When 'x' is an energy in meV, the parameter p should be T/11.6045
%
% Reference: http://en.wikipedia.org/wiki/Bose-Einstein_statistics
%
% input:  p: Bose model parameters (double)
%            p = [ Kb.T/hbar   in 'x' units ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=bose(300/11.605, -10:10); or y=plot(bose);
%
% Version: $Revision$
% See also iData, iFunc/fits, iFunc/plot

y.Name       = [ 'Bose (1D) [' mfilename ']' ];
y.Description='Bose-Einstein distribution fitting function. Ref: http://en.wikipedia.org/wiki/Bose%E2%80%93Einstein_statistics';
y.Parameters = {'Tau h/2pi/kT'};
y.Expression = @(p,x) 1 ./ (exp(abs(x)/p(1)) - 1) + (x>0);
y.Dimension  = 1;   
y.Guess      = @(x,signal) abs(log(1./mean(signal(:))+1)/mean(abs(x(:))));

y = iFunc(y);

if length(varargin)
  y = y(varargin{:});
end
