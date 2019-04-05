function y=sf_square_well(varargin)
% y = sf_square_well(p, x, [y]) : structure factor of particles interacting with a square well potential [Sharma] S(q)
%
%   iFunc/sf_square_well structure factor S(q) of particles interacting with a 
%       square well potential [Sharma]
%     The 'x' wave-vector/momentum axis is usually in nm-1 or Angs-1.
%     The parameter 'R' is given in inverse unit of the axis (that is nm or Angs)
%       and corresponds with the typical distance between scattering objects.
%     The parameter 'rho' is a reduced density. A value of rho close to 1 is for 
%       a cristalline material. A value of 0 corresponds with a perfect liquid/gas S(q)=1.
%     Typical values for parameters are 
%       R=10 Angs, rho=0.02, epsilon=1.5*T, Delta=1.2*R
%     The model returns the S(q) structure factor.
%
%     Ref: Sharma, R. V.; Sharma, K. C. Physica, 89A, 213. (1977).
%          Extracted from sasfit/sasfit_sq/sasfit_sq_SquareWell2.c
%     I. Bressler, et al, Journal of Applied Crystallography, 2015, 48 (5), 1587-1598
%
% input:  p: square well potential model parameters (double)
%            p = [ R=Hard_Sphere_Radius rho=Volume_Fraction epsilon=Well_Depth_K Delta=WellWidth T=Temperature_K]
%          or 'guess'
%         x: wave-vector/momentum axis (double, e.g. nm-1 or Angs-1)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value (intensity S(q))
% ex:     y=sf_square_well([50 0.05 15 10 10], -10:10); or plot(sf_square_well,[10 0.02 15 10 10], 0:0.01:1)
%
% Version: $Date$ $Version$ $Author$
% See also iFunc, iFunc/fits, iFunc/plot
% 

y.Name      = [ 'Square well S(q) (1D) [' mfilename ']' ];
y.Description='Square well potential scattering structure factor [Sharma]';
y.Parameters={'R hard sphere radius [1/x]', ...
              'rho hard sphere volume fraction', ...
              'epsilon square well depth [K]', ...
              'Delta square well width [1/x]', ...
              'T Temperature [K]'};
y.Dimension =1;
y.Expression= { ...
  'RHS=p(1); eta=p(2); epsi_div_kB=p(3); Delta=p(4); T=p(5);q=abs(x);' ...
  'kb = 1.3806505E-23;' ...
	'epsilon = epsi_div_kB*kb;' ...
	'sigma = 2*RHS;' ...
	'sigma3 = sigma*sigma*sigma;' ...
	'rho = eta*6/pi/sigma3;' ...
	'kappa = sigma*q;' ...
	'kappa3=kappa.^3;' ...
	'signal = 1.-4*pi*rho*sigma3*(sin(kappa)-kappa.*cos(kappa))./kappa3;' ...
	'signal = signal+4*pi*rho*sigma*sigma*(exp(epsilon/(kb*T))-1.)*sin(kappa)./kappa*Delta;'
};
% use ifthenelse anonymous function
% <https://blogs.mathworks.com/loren/2013/01/10/introduction-to-functional-programming-with-anonymous-functions-part-1/>
% iif( cond1, exec1, cond2, exec2, ...)
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
y.Guess     = @(x,signal) iif(...
     ~isempty(signal)&&numel(signal)==numel(x), @() [ pi/sum(signal(:).*x(:))*sum(signal(:)) max(max(signal(:)-1),0.01) ...
        15 1.2*pi/sum(signal(:).*x(:))*sum(signal(:)) 10 ], ...
     true,             @() [50 0.05 15 10 10]);
y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

