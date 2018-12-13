function signal=sab_diffusion(varargin)
% model = sab_diffusion(p, alpha ,beta, {signal}) : Free-gas/recoil dispersion(Q)
%
%   iFunc/sab_diffusion: a 2D S(alpha,beta) with a diffusion dispersion
%     based on the Egelstaff-Schofield "effective width model".
%     This is a pure incoherent scattering law (no structure).
%
%   The input parameter 'wt' is the translational weight, which quantifies the 
%   fraction of the scattering originating from recoil/diffusion, e.g. wt=Mdiff/M
%   where Mdiff is the mass that diffuses. A pure diffusion is obtained with wt=1.
%   The input parameter 'c' quantifies the diffusion coefficient, c=MD/wt/h. 
%   A standard recoil/free-gas is obtained for large 'c' values (e.g. 4).
%   An incoherent solid type diffusion is obtained for small 'c' values (e.g. 0.01).
%   The model satisfies the detailed balance.
%
%   The dispersion has the form:
%      S(a,b) = 2*c*wt/pi*a * exp(2*c^2*a-b/2) 
%               * sqrt( (c^2+1/4)./(b^2+4*c^2*wt^2*a^2) )
%               * K1( sqrt( (c^2+0.25)*(b^2+4*c^2*wt^2*a^2) ) )
%
%   where we commonly define (h stands for hbar):
%     a    =  h2q2/2MkT = (Ei+Ef-2*mu*sqrt(Ei*Ef))/AkT    unit-less moment
%     b    = -hw/kT     = (Ef-Ei)/kT                      unit-less energy
%     A    = M/m                
%     mu   = cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)
%     m    = mass of the neutron
%     M    = mass of the scattering target
%     Mdiff= mass of the diffusive scattering target, Mdiff <= M
%     wt   = translational weight, wt=Mdiff/M e.g. in [0-1]
%     c    = unitless diffusion coefficient, c=MD/wt/h, e.g. [0-4]
%
% conventions:
% w = omega = Ei-Ef = energy lost by the neutron [meV]
%     omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%     omega < 0, neutron gains energy, anti-Stokes
%
%   You can build a diffusive model for a given translational weight and diffusion
%   coeeficient:
%      sab = sab_diffusion([ wt c ])
%
%   You can of course tune other parameters once the model object has been created.
%
%   Evaluate the model with syntax:
%     sab(p, alpha, beta)
%
% input:  p: sab_diffusion model parameters (double)
%             p(1)= wt        Translational weight wt=Mdiff/M [1]
%             p(2)= c         Unitless diffusion constant, c=MD/h [typically 0-4]
%             p(3)= Amplitude 
%         alpha:  axis along unitless wavevector/momentum (row,double)
%         beta:   axis along unitless energy (column,double)
% output: signal: model value [iFunc_Sab]
%
%
% Example:
%   s=sab_diffusion;
%   plot(log10(iData(s, [], 0:.1:20, -50:50)))  % alpha=0:20, beta=-50:50
%
% Reference: 
%  Egelstaff and Schofield, Nuc. Sci. Eng. 12 (1962) 260 <https://doi.org/10.13182/NSE62-A26066>
%  J.I. Marquez-Damian et al, Ann. Nuc. En. 92 (2016) 107 <http://dx.doi.org/10.1016/j.anucene.2016.01.036>
%  M.Mattes and J.Keinert, IAEA INDC (NDS)-0470 (2005) https://www-nds.iaea.org/publications/indc/indc-nds-0470/
%  R.E.McFarlane, LA-12639-MS (ENDF 356) (March 1994) https://t2.lanl.gov/nis/publications/thermal.pdf
%  
%
% Version: $Date$
% See also iData, iFunc
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

signal.Name           = [ 'sab_diffusion dispersion(Q) free-gas dispersion [' mfilename ']' ];
signal.Description    = 'A 2D S(alpha,beta) free-gas/translational dispersion.';

signal.Parameters     = {  ...
  'wt             Translational weight; wt=Mdiff/M [1]' ...
  'c              Diffusion constant; c=MD/wt/h, unitless [e.g. 0-4]' ...
  'Amplitude' };
  
signal.Dimension      = 2;         % dimensionality of input space (axes) and result
signal.Guess          = [ 1 1 1 ];

signal.Expression     = { ...
 'a = x; b = y; wt = p(1); c = p(2);' ...
 'if isvector(a) && isvector(b) && numel(a) ~= numel(b), [a,b] = meshgrid(a,b); end' ...
 'signal = p(3)*2*c*wt/pi*a .* exp(2*c^2*a-b/2) .* sqrt( (c^2+0.25)./(b.^2+4*c^2*wt^2.*a.^2) ) .* besselk(1, sqrt( (c^2+0.25)*(b.^2+4*c^2*wt^2*a.^2) ) );' ...
 };

signal= iFunc(signal);
signal= iFunc_Sab(signal); % overload Sab flavour

if nargin == 1 && isnumeric(varargin)
  p = varargin{1};
  if numel(p) >= 1, signal.ParameterValues(1) = p(1); end
  if numel(p) >= 2, signal.ParameterValues(2) = p(2); end
end
