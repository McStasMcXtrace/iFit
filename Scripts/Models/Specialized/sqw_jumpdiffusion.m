function signal=sqw_jumpdiffusion(varargin)
% model = sqw_jumpdiffusion(p, q ,w, {signal}) : jump/Fick-law diffusion dispersion(Q) Sqw2D
%
%   iFunc/sqw_jumpdiffusion: a 2D S(q,w) with a jump diffusion dispersion
%     based on the Egelstaff model. It can also model a simple Fick's law.
%     This is a classical pure incoherent Lorentzian scattering law (no structure).
%
%  Model and parameters:
%  ---------------------
%
%   The dispersion has the form: (Egelstaff book Eq 11.13 and 11.16, p 222)
%      S(q,w) = f(q)/(w^2+f(q)^2)
%      f(q)   = w0 q^2 l0^2/(1+q^2 l0^2)
%
%   where we commonly define:
%     w0   = Jump diffusion characteristic energy width, MD/2kT e.g. few [meV]
%     l0   = Jump diffusion length e.g. few [Angs].
%
%   The characteristic energy for a jump step is w0, usually around 
%   few meV in liquids, which inverse time t0 characterises the residence time 
%   step between jumps, t0 ~ 1-4 ps. The half width of the Lorentizian
%   follows a law similar to Dq^2.
%
%   The mean free path is l0 is around 0.1-5 Angs. 
%
%   When q l0 is small, f(q) -> w0 q^2 l0^2 which behaves as a free-diffusion
%   Fick's law with an equivalent diffusion constant D=w0.l0^2 usally around 
%   D=1-10 E-9 m^2/s in liquids.
%
%   You can build a jump diffusion model for a given translational weight and 
%   diffusion coefficient:
%      sqw = sqw_jumpdiffusion([ w0 l0 ])
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
%    sqw = sqw_diffusion .* bose;
%  where the Temperature is then given in [x units]. If 'x' is an energy in [meV]
%  then the Temperature parameter is T[K]/11.6045
%
%  Energy conventions:
%   w = omega = Ei-Ef = energy lost by the neutron [meV]
%       omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%       omega < 0, neutron gains energy, anti-Stokes
%
% Usage:
% ------
% s = sqw_jumpdiffusion; 
% value = s(p,q,w); or value = iData(s,p,q,w)
%
% input:  p: sqw_jumpdiffusion model parameters (double)
%             p(1)= Amplitude 
%             p(2)= w0             Jump diffusion characteristic energy width, e.g. few [meV]
%             p(2)= l0             Jump diffusion length [Angs]
%         q:  axis along wavevector/momentum (row,double) [Angs-1]
%         w:  axis along energy (column,double) [meV]
% output: signal: model value [iFunc_Sqw2D]
%
% Example:
%   s=sqw_jumpdiffusion;
%   plot(log10(iData(s, [], 0:.1:20, -50:50)))  % q=0:20 Angs-1, w=-50:50 meV
%
% Reference: 
%  P.A.Egelstaff, An introduction to the liquid state, 2nd ed., Oxford (2002)
%  Egelstaff and Schofield, Nuc. Sci. Eng. 12 (1962) 260 <https://doi.org/10.13182/NSE62-A26066>
%  J.I. Marquez-Damian et al, Ann. Nuc. En. 92 (2016) 107 <http://dx.doi.org/10.1016/j.anucene.2016.01.036>
%
% Version: $Date$
% See also iData, iFunc, sqw_recoil, sqw_diffusion
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

signal.Name           = [ 'sqw_jumpdiffusion jump/Fick-law diffusion dispersion [' mfilename ']' ];
signal.Description    = 'A 2D S(q,w) jump diffusion dispersion.';

signal.Parameters     = {  ...
  'Amplitude' ...
  'w0             Jump diffusion characteristic energy width [meV]' ...
  'l0             Jump diffusion length [Angs]' ...
   };
  
signal.Dimension      = 2;         % dimensionality of input space (axes) and result
signal.Guess          = [ 1 1 1 ];
signal.UserData.classical     = true;
signal.UserData.DebyeWaller   = false;

% Expression of S(q,w) is found in Egelstaff, Intro to Liquids, Eq (11.13)+(11.16), p 227.
signal.Expression     = { ...
 'q = x; w = y; w0 = p(2); l0 = p(3);' ...
 'if isvector(q) && isvector(w) && numel(q) ~= numel(w), [q,w] = meshgrid(q,w); end' ...
 'fq = w0*l0^2*q.^2./(1+l0^2*q.^2); % half width of Lorentzian ~ DQ^2 in [meV]' ...
 'signal = p(1)/pi*fq./(w.^2+fq.^2);' ...
 };

signal= iFunc(signal);
signal= iFunc_Sqw2D(signal); % overload Sqw flavour

if nargin == 1 && isnumeric(varargin{1})
  p = varargin{1};
  if numel(p) >= 1, signal.ParameterValues(2) = p(1); end
  if numel(p) >= 2, signal.ParameterValues(3) = p(2); end
end
