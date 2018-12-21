function signal=sqw_diffusion(varargin)
% model = sqw_diffusion(p, q ,w, {signal}) : Brownian diffusion dispersion(Q)
%
%   iFunc/sqw_diffusion: a 2D S(q,w) with a diffusion dispersion
%     based on the Egelstaff-Schofield Langevin equation for Brownian motion model.
%     This is a classical pure incoherent Gaussian scattering law (no structure).
%
%   The model is classical, e.g. symmetric in energy, and does not satisfy the
%   detailed balance.
%
%   The characteristic energy for a translational diffusion step is w0, usually  
%   around few meV in liquids, corresponds to a friction, which inverse time t0
%   characterises the diffusion step, t0 ~ 1-4 ps. Usually, one can write
%   w0 = 1/t0 = MD/2kT.
%
%   The mean free path is computed as l0 = sqrt(6*t0*D) and is around 0.1-0.5 nm.
%   The diffusion constant D is usally around D=1-10 E-9 m^2/s in liquids.
%
%   The dispersion has the form: (Egelstaff book Eq 11.32, p 227.)
%
%      S(q,w) = exp(Dq^2/w0) Dq^2/w0/sqrt(w^2+(Dq^2)^2) K1(sqrt(w^2+(Dq^2)^2)/w0)
%
%   where K1 is a modified Bessel function of the second kind.
%
%   where we commonly define:
%     w0   = Diffusion characteristic energy, MD/2kT e.g. few [meV]
%     D    = Diffusion coefficient e.g. few 1E-9 [m^2/s]
%
% conventions:
% w = omega = Ei-Ef = energy lost by the neutron [meV]
%     omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%     omega < 0, neutron gains energy, anti-Stokes
%
%   You can build a diffusive model for a given diffusion energy and
%   coefficient:
%      sqw = sqw_diffusion([ w0 D ])
%
%   You can of course tune other parameters once the model object has been created.
%
%   Evaluate the model with syntax:
%     sqw(p, q, w)
%
% input:  p: sqw_diffusion model parameters (double)
%             p(1)= Amplitude 
%             p(2)= w0             Diffusion characteristic energy [meV], e.g. 1
%             p(2)= D              Diffusion constant [m^2/s] e.g. 1E-9
%         q:  axis along wavevector/momentum (row,double) [Angs-1]
%         w:  axis along energy (column,double) [meV]
% output: signal: model value [iFunc_Sqw2D]
%
%
% Example:
%   s=sqw_diffusion;
%   plot(log10(iData(s, [], 0:.1:20, -50:50)))  % q=0:20 Angs-1, w=-50:50 meV
%
% Reference: 
%  P.A.Egelstaff, An introduction to the liquid state, 2nd ed., Oxford (2002)
%  Egelstaff and Schofield, Nuc. Sci. Eng. 12 (1962) 260 <https://doi.org/10.13182/NSE62-A26066>
%  J.I. Marquez-Damian et al, Ann. Nuc. En. 92 (2016) 107 <http://dx.doi.org/10.1016/j.anucene.2016.01.036>
%  M.Mattes and J.Keinert, IAEA INDC (NDS)-0470 (2005) https://www-nds.iaea.org/publications/indc/indc-nds-0470/
%  R.E.McFarlane, LA-12639-MS (ENDF 356) (March 1994) https://t2.lanl.gov/nis/publications/thermal.pdf
%
% Version: $Date$
% See also iData, iFunc, sqw_recoil, sab_diffusion, sab_recoil
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

signal.Name           = [ 'sqw_diffusion Brownian diffusion dispersion [' mfilename ']' ];
signal.Description    = 'A 2D S(q,w) Brownian diffusion dispersion.';

signal.Parameters     = {  ...
  'Amplitude' ...
  'w0             Diffusion characteristic energy width [meV]' ...
  'D              Diffusion constant [m^2/s]' ...
   };
  
signal.Dimension      = 2;         % dimensionality of input space (axes) and result
signal.Guess          = [ 1 1 1E-9 ];

% Expression of S(q,w) is found in Egelstaff, Intro to Liquids, Eq (11.32), p 227.

% w0 = kT/MD
signal.Expression     = { ...
 'q = x; w = y; w0 = p(2); D = p(3);' ...
 'if isvector(q) && isvector(w) && numel(q) ~= numel(w), [q,w] = meshgrid(q,w); end' ...
 'dq2    = D*q.^2*1E20/241.8E9; % in meV, with q in Angs-1' ...
 'wdq2   = sqrt(w.^2+dq2.^2);   % in meV' ...
 'signal = p(1)*exp(dq2/w0) .* dq2/w0./wdq2 .* besselk(1, wdq2/w0);' ...
 };

signal= iFunc(signal);
signal= iFunc_Sqw2D(signal); % overload Sqw flavour

if nargin == 1 && isnumeric(varargin{1})
  p = varargin{1};
  if numel(p) >= 1, signal.ParameterValues(2) = p(1); end
  if numel(p) >= 2, signal.ParameterValues(3) = p(2); end
end
