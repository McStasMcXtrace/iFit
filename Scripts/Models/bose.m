function y=bose(varargin)
% y = bose(p, x, [y]) : Bose
%
%   iFunc/bose Bose fitting function
%     n(x) = 1 ./ (exp(x/p(1)) - 1);
%     y = n(x)+1 for x > 0
%     y = n(x)   for x < 0
%
% The Bose factor has the sign of x.
%   (n+1) converges to 0 for x -> -Inf, and to 1 for x-> +Inf. It diverges at x=0.
%   n(x)+n(-x) = -1
%         
% When 'x' is an energy in meV, the parameter p should be T/11.6045
%
% Energy conventions:
%   x = omega = Ei-Ef = energy lost by the neutron [meV]
%       omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%       omega < 0, neutron gains energy, anti-Stokes
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
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot
% (c) E.Farhi, ILL. License: EUPL.

y.Name       = [ 'Bose (1D) [' mfilename ']' ];
y.Description='Bose-Einstein distribution fitting function. Ref: http://en.wikipedia.org/wiki/Bose%E2%80%93Einstein_statistics';
y.Parameters = {'Temperature [x unit, 2pi*kT/h]'};
% the Bose factor is negative for w<0, positive for w>0
% (n+1) converges to 0 for w -> -Inf, and to 1 for w-> +Inf. It diverges at w=0
y.Expression = { ...
  'signal = 1 ./ (exp(x/p(1)) - 1);' ...
  'signal(find(x==0)) = 1;' ...
};
y.Dimension  = 1;   
% use ifthenelse anonymous function
% <https://blogs.mathworks.com/loren/2013/01/10/introduction-to-functional-programming-with-anonymous-functions-part-1/>
% iif( cond1, exec1, cond2, exec2, ...)
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
y.Guess     = @(x,signal) iif(...
  ~isempty(signal), @() abs(log(1./mean(signal(:))+1)/mean(abs(x(:)))), ...
  true            , @() 25);

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end
