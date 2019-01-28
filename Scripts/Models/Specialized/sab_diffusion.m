function signal=sab_diffusion(varargin)
% model = sab_diffusion(p, alpha ,beta, {signal}) : Brownian diffusion dispersion(Q) Sab
%
%   iFunc/sab_diffusion: a 2D S(alpha,beta) with a diffusion dispersion
%     based on the Egelstaff-Schofield Langevin equation for Brownian motion model.
%     This is a pure incoherent Gaussian scattering law (no structure).
%     This model is equivalent to the NJOY/LEAPR "diffusion effective width model".
%
%   The input parameter 'wt' is the translational weight, which quantifies the 
%   fraction of the scattering originating from recoil/diffusion, e.g. wt=Mdiff/M
%   where Mdiff is the mass that diffuses. A pure diffusion is obtained with wt=1.
%   The input parameter 'c' quantifies the diffusion coefficient, c=MD/wt/h. 
%   A standard recoil/free-gas is obtained for large 'c' values (e.g. 4).
%   An incoherent solid type diffusion is obtained for small 'c' values (e.g. 0.01).
%   The model satisfies the detailed balance.
%
%   The characteristic time for a diffusion step is t0=MD/kT, usually around 
%   t0 = 1-4 1E-12 s in liquids.
%   The diffusion constant D is usally around D= 1-10  1E-3 mm2/s in liquids.
%
%   The dispersion has the form (Egelstaff and Schofield NSE 1962, Eq 4.8):
%      S(a,b) = 2*c*wt/pi*a * exp(2*c^2*a-b/2) 
%               * sqrt( (c^2+1/4)./(b^2+4*c^2*wt^2*a^2) )
%               * K1( sqrt( (c^2+0.25)*(b^2+4*c^2*wt^2*a^2) ) )
%
%   where we commonly define (h stands for hbar):
%     a    =  h2q2/2MkT = (Ei+Ef-2*mu*sqrt(Ei*Ef))/AkT    unit-less moment alpha
%     b    = -hw/kT     = (Ef-Ei)/kT                      unit-less energy beta
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
%   coefficient:
%      sab = sab_diffusion([ wt c ])
%
%   You can of course tune other parameters once the model object has been created.
%
%   Evaluate the model with syntax:
%     sab(p, alpha, beta)
%
% input:  p: sab_diffusion model parameters (double)
%             p(1)= Amplitude 
%             p(2)= wt        Translational weight wt=Mdiff/M [1]
%             p(3)= c         Unitless diffusion constant, c=MD/h [typically 0-4]
%         alpha:  axis along unitless wavevector/momentum (row,double)
%         beta:   axis along unitless energy (column,double)
% output: signal: model value [iFunc_Sab]
%
% Example:
%   s=sab_diffusion;
%   plot(log10(iData(s, [], 0:.1:20, -50:50)))  % alpha=0:20, beta=-50:50
%
% Reference: 
%  P.A.Egelstaff, An introduction to the liquid state, 2nd ed., Oxford (2002)
%  Egelstaff and Schofield, Nuc. Sci. Eng. 12 (1962) 260 <https://doi.org/10.13182/NSE62-A26066>
%  J.I. Marquez-Damian et al, Ann. Nuc. En. 92 (2016) 107 <http://dx.doi.org/10.1016/j.anucene.2016.01.036>
%  M.Mattes and J.Keinert, IAEA INDC (NDS)-0470 (2005) https://www-nds.iaea.org/publications/indc/indc-nds-0470/
%  R.E.McFarlane, LA-12639-MS (ENDF 356) (March 1994) https://t2.lanl.gov/nis/publications/thermal.pdf
%  
% Version: $Date$
% See also iData, iFunc, sqw_diffusion, sab_recoil, sqw_recoil
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

signal.Name           = [ 'sab_diffusion Brownian diffusion dispersion [' mfilename ']' ];
signal.Description    = 'A 2D S(alpha,beta) Brownian diffusion dispersion.';

signal.Parameters     = {  ...
  'Amplitude' ...
  'wt             Translational weight; wt=Mdiff/M [1]' ...
  'c              Diffusion constant; c=MD/wt/h, unitless [e.g. 0-4]' ...
   };
  
signal.Dimension      = 2;         % dimensionality of input space (axes) and result
signal.Guess          = [ 1 1 1 ];
signal.UserData.classical     = false;
signal.UserData.DebyeWaller   = true;

% Egelstaff and Schofield NSE (1962) , Eq (4.8)
signal.Expression     = { ...
 'a = x; b = y; wt = p(2); c = p(3);' ...
 'signal = p(1)*2*c*wt/pi*a .* exp(2*c^2*a-b/2) .* sqrt( (c^2+0.25)./(b.^2+4*c^2*wt^2.*a.^2) ) .* besselk(1, sqrt( (c^2+0.25)*(b.^2+4*c^2*wt^2*a.^2) ) );' ...
 };

signal= iFunc(signal);
signal= iFunc_Sab(signal); % overload Sab flavour

if nargin == 1 && isnumeric(varargin{1})
  p = varargin{1};
  if numel(p) >= 1, signal.ParameterValues(2) = p(1); end
  if numel(p) >= 2, signal.ParameterValues(3) = p(2); end
end
