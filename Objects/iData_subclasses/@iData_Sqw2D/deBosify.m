function s = deBosify(s, T, type)
% iData_Sqw2D: deBosify: remove Bose factor (detailed balance) from an 'experimental/quantum' data set.
%   The initial data set should be 'quantum/experimental' and satisfy the detailed balance.
%   The resulting data set obeys S*=S(q,w) = S(q,-w), i.e. is 'classical'.
%   It suppresses the temperature effect (population).
%
% conventions:
% omega = Ei-Ef = energy lost by the neutron, given in [meV]
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
%
%    S(q,-w) = exp(-hw/kT) S(q,w)
%    S(q,w)  = exp( hw/kT) S(q,-w)
%    S(q,w)  = Q(w) S*(q,w) with S*=classical limit and Q(w) defined below.
% for omega > 0, S(q,w) > S(q,-w)
%               
% The semi-classical correction, Q, aka 'quantum' correction factor, 
% can be selected from the optional   'type' argument:
%    Q = exp(hw_kT/2)                 'Schofield' or 'Boltzmann'
%    Q = hw_kT./(1-exp(-hw_kT))       'harmonic'  or 'Bader'
%    Q = 2./(1+exp(-hw_kT))           'standard'  or 'Frommhold' (default)
%
% The 'Boltzmann' correction leads to a divergence of the S(q,w) for e.g. w above 
% few 100 meV. The 'harmonic' correction provides a reasonable correction but does
% not fully avoid the divergence at large energies.
%
%  Bose factor: n(w) = 1./(exp(w*11.605/T) -1) ~ exp(-w*11.605/T)
%               w in [meV], T in [K]
%
% References:
%  B. Hehr, http://www.lib.ncsu.edu/resolver/1840.16/7422 PhD manuscript (2010).
%  S. A. Egorov, K. F. Everitt and J. L. Skinner. J. Phys. Chem., 103, 9494 (1999).
%  P. Schofield. Phys. Rev. Lett., 4, 239 (1960).
%  J. S. Bader and B. J. Berne. J. Chem. Phys., 100, 8359 (1994).
%  T. D. Hone and G. A. Voth. J. Chem. Phys., 121, 6412 (2004).
%  L. Frommhold. Collision-induced absorption in gases, 1 st ed., Cambridge
%    Monographs on Atomic, Molecular, and Chemical Physics, Vol. 2,
%    Cambridge Univ. Press: London (1993).
%
% syntax:
%   sqw = deBosify(sqw_T)
%
% input:
%   sqw_T: Sqw data set (non classical, including T Bose factor e.g from experiment)
%           e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
%   T:     when given, Temperature to use for Bose. When not given, the Temperature
%           is searched in the object.
%   type: 'Schofield' or 'harmonic' or 'standard' (default)
%
% output:
%   sqw: classical/symmetric Sqw data set (iData_Sqw2D).
%
% Example: s=iData_Sqw2D('SQW_coh_lGe.nc'); sb=Bosify(s, 1235); s0=deBosify(sb);
%
% See also: iData_Sqw2D/Bosify, iData_Sqw2D/symmetrize, 
%           iData_Sqw2D/dynamic_range, iData_Sqw2D/scattering_cross_section
% (c) E.Farhi, ILL. License: EUPL.
  
  if nargin < 2, T = []; end
  if nargin < 3, type=''; end
  if isempty(type), type='standard'; end
  
  s = Bosify(s, T, [ type ' debosify' ]);
  

  
