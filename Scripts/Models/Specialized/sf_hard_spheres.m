function y=sf_hard_spheres(varargin)
% y = sf_hard_spheres(p, x, [y]) : Hard Sphere structure factor [Percus-Yevick] S(q)
%
%   iFunc/sf_hard_spheres Hard Sphere structure factor S(q), suited for simple liquids.
%     It corresponds to hard spheres in a Lennard-Jones potential.
%     The 'x' wave-vector/momentum axis is usually in nm-1 or Angs-1.
%     The parameter 'R' is given in inverse unit of the axis (that is nm or Angs)
%       and corresponds with the typical distance between scattering objects.
%     The parameter 'rho' is a reduced density. A value of rho close to 1 is for 
%       a cristalline material. A value of 0 corresponds with a perfect liquid/gas S(q)=1.
%
%     Typical values for parameters are R=3-50 Angs, rho=0.4.
%     The Hard Sphere model corresponds with the Sticky Hard Sphere model with large tau.
%     The model returns the S(q) structure factor.
%
%     Ref: J. K. Percus and G. J. Yevick., Phys. Rev., 110(1):1–13, 1958.
%          A. Vrij., J. Chem. Phys., 71(8):3267–3270, 1979.
%          Extracted from sasfit/sasfit_sq/sasfit_sq_HardSphere.c
%     I. Bressler, et al, Journal of Applied Crystallography, 2015, 48 (5), 1587-1598
%
% input:  p: hard sphere model parameters (double)
%            p = [ R=Hard_Sphere_Radius rho=Volume_Fraction ]
%          or 'guess'
%         x: wave-vector/momentum axis (double, e.g. nm-1 or Angs-1)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value (intensity)
% ex:     y=sf_hard_spheres([4 0.4], 0:0.01:1); or plot(sf_hard_spheres,[10 0.1],0:0.01:1)
%
% Version: $Date$
% See also iFunc, iFunc/fits, iFunc/plot
% (c) E.Farhi, ILL. License: EUPL.

y.Name      = [ 'Hard_Sphere S(q) (1D) [' mfilename ']' ];
y.Description='Hard Sphere scattering structure factor [Percus-Yevick]';
y.Parameters={'R hard sphere radius [1/x]', ...
              'rho hard sphere volume fraction'};
y.Expression= { ...
  'if exist(''y'')' ...
  'if isvector(x) && isvector(y) && ~isempty(y) && numel(x) ~= numel(y), [x,y] = meshgrid(x,y); end' ...
  'end' ...
  'A=2.0*abs(p(1)*x); fp=max(0, min(p(2), 1)); '...
  'if (p(2) <= 0.0) signal=ones(size(x)); else ' ...
  'alpha = power(1.0+2.0*fp,2.0)/power(1.0-fp,4.0); ' ...
  'beta  = -6.0*fp*power(1.0+fp/2.0,2.0)/power(1.0-fp,4.0); gamma = fp*alpha/2.0; ' ...
  'signal = alpha*(sin(A)-A.*cos(A))./power(A,2.0); ' ...
  'signal = signal + beta*(2.0*A.*sin(A)+(2.0-power(A,2.0)) .* cos(A)-2.0)./power(A,3.0); ' ...
  'signal = signal + gamma * (-power(A,4.0) .* cos(A) + 4.0*((3.0*power(A,2.0)-6.0).*cos(A)+(power(A,3.0)-6.0*A).*sin(A)+6.0))./power(A,5.0); ' ...
  'signal = 1.0./(1.0+24.0*fp*signal./A); signal(isnan(signal))=0; end'
};
% use ifthenelse anonymous function
% <https://blogs.mathworks.com/loren/2013/01/10/introduction-to-functional-programming-with-anonymous-functions-part-1/>
% iif( cond1, exec1, cond2, exec2, ...)
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
y.Guess     = @(x,signal) iif(...
  ~isempty(signal), @() [ pi/sum(signal(:).*x(:))*sum(signal(:)) max(max(signal(:)-1),0.01) ], ...
  true            , @() [ 4 0.4 ]);
y.Dimension = 1;
y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

