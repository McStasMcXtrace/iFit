function signal=sqw_diffusion_jump(varargin)
% model = sqw_diffusion_jump(p, q ,w, {signal}) : jump/Fick-law diffusion dispersion(Q) Sqw2D
%
%   iFunc/sqw_diffusion_jump: a 2D S(q,w) with a jump diffusion dispersion
%     based on the Egelstaff model. It can also model a simple Fick's law.
%     It is suited for quasi-elastic neutron scattering data (QENS).
%     This is a classical pure incoherent Lorentzian scattering law (no structure).
%
%  Model and parameters:
%  ---------------------
%
%   The dispersion has the form: (Egelstaff book Eq 11.13 and 11.16, p 222)
%      S(q,w) = f(q)/(w^2+f(q)^2)
%   where
%      f(q)   = Dq^2 /(1+t0.Dq^2)     width of the Lorentzian
%
%   where we commonly define:
%     D    = Diffusion coefficient e.g. few 1E-9 [m^2/s]
%     t0   = Residence time step between jumps, e.g. 1-4 ps [s]
%
%   The diffusion constant D is usually around D=1-10 E-9 [m^2/s] in liquids. Its 
%   value in [meV.Angs^2] is D*4.1356e+08. The residence time t0 is usually in 0-4 ps.
%
%   Fixing t0=0 corresponds with a pure random walk translation (Brownian/Ficks-law
%   self-diffusion). Liquid particles then collide randomly (directional memory loss).
%   The Brownian mean free path l0=sqrt(6*t0*D) is usually around 0.1-5 Angs.
%
%   You can build a jump diffusion model for a given translational weight and 
%   diffusion coefficient:
%      sqw = sqw_diffusion_jump([ D t0 ])
%
%   You can of course tune other parameters once the model object has been created.
%
%   Evaluate the model with syntax:
%     sqw(p, q, w)
%
%  Additional remarks:
%  -------------------
%
%  The model is classical, e.g. symmetric in energy, and does NOT satisfy the
%  detailed balance.
%
%  To get the 'true' quantum S(q,w) definition, use e.g.
%    sqw = Bosify(sqw_diffusion_jump);
%
%  To add a Debye-Waller factor (thermal motions around the equilibrium), use e.g.
%    sqw = DebyeWaller(sqw);
%
%  Energy conventions:
%   w = omega = Ei-Ef = energy lost by the neutron [meV]
%       omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%       omega < 0, neutron gains energy, anti-Stokes
%
%  This model is equivalent to the QENS_JumpDiffusionTranslation model in LAMP.
%  and QENS_RandomWalkTranslation when t0=0 [http://lamp.mccode.org].
%
% Usage:
% ------
% s = sqw_diffusion_jump; 
% value = s(p,q,w); or value = iData(s,p,q,w)
%
% input:  p: sqw_diffusion_jump model parameters (double)
%             p(1)= Amplitude 
%             p(2)= D         Diffusion coefficient e.g. few 1E-9 [m^2/s]
%             p(3)= t0        Residence time step between jumps, e.g. 1-4 ps [s]
%         q:  axis along wavevector/momentum (row,double) [Angs-1]
%         w:  axis along energy (column,double) [meV]
% output: signal: model value [iFunc_Sqw2D]
%
% Example:
%   s=sqw_diffusion_jump;
%   plot(log10(iData(s, [], 0:.1:20, -50:50)))  % q=0:20 Angs-1, w=-50:50 meV
%
% Reference: 
%  P.A.Egelstaff, An introduction to the liquid state, 2nd ed., Oxford (2002)
%  Egelstaff and Schofield, Nuc. Sci. Eng. 12 (1962) 260 <https://doi.org/10.13182/NSE62-A26066>
%  J.I. Marquez-Damian et al, Ann. Nuc. En. 92 (2016) 107 <http://dx.doi.org/10.1016/j.anucene.2016.01.036>
%  M. Bée: Quasielastic Neutron Scattering: Principles and Applications in Solid State Chemistry, Biology and Materials Science, Adam-Hilger, Bristol, 1988; chapter 5.1.
%  R. Hempelman: Quasielastic Neutron Scattering and Solid State Diffusion, Clarendon Press, Oxford, 2000; chapter 5.2
%
% Version: $Date$
% See also iData, iFunc, sqw_recoil, sqw_diffusion
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

signal.Name           = [ 'sqw_diffusion_jump jump/Fick-law diffusion dispersion [' mfilename ']' ];
signal.Description    = 'A 2D S(q,w) jump diffusion dispersion.';

signal.Parameters     = {  ...
  'Amplitude' ...
  'D              Diffusion coefficient e.g. few 1E-9 [m^2/s]' ...
  't0             Residence time step between jumps, e.g. 1-4 ps [s]' ...
   };
  
signal.Dimension      = 2;         % dimensionality of input space (axes) and result
signal.Guess          = [ 1 1e-9 1e-12 ];
signal.UserData.classical     = true;
signal.UserData.DebyeWaller   = false;

% Expression of S(q,w) is found in Egelstaff, Intro to Liquids, Eq (11.13)+(11.16), p 227.
signal.Expression     = { ...
 'q = x; w = y; D = p(2); t0 = p(3)*241.8E9; % in meV-1' ...
 'dq2= D*q.^2*1E20/241.8E9; % in meV, with q in Angs-1' ...
 '% half width of Lorentzian ~ DQ^2 in [meV]' ...
 'fq=dq2;' ...
 'if t0>0, fq = fq./(1+t0*dq2); end' ...
 'signal = p(1)/pi*fq./(w.^2+fq.^2);' ...
 };

signal= iFunc(signal);
signal= iFunc_Sqw2D(signal); % overload Sqw flavour

if nargin == 1 && isnumeric(varargin{1})
  p = varargin{1};
  if numel(p) >= 1, signal.ParameterValues(2) = p(1); end
  if numel(p) >= 2, signal.ParameterValues(3) = p(2); end
end
