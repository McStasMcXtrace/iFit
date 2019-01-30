function self=deBosify(self, type)
% iFunc_Sqw2D: deBosify: remove the 'Bose' factor (detailed balance) to a classical S(q,w) model.
%   The initial model should be 'quantum/experimental' and satisfy the detailed balance.
%   The resulting model obeys S*=S(q,w) = S(q,-w), i.e. is 'classical'.
%   It suppresses most of the temperature effect (population).
%
% conventions:
% omega = w = Ei-Ef = energy lost by the neutron, given in [meV]
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
%
%    S(q,-w) = exp(-hw/kT) S(q,w)
%    S(q,w)  = exp( hw/kT) S(q,-w)
%    S(q,w)  = Q(w) S*(q,w) with S*=classical limit and Q(w) defined below.
%    for omega > 0, S(q,w) > S(q,-w)
%         
% The semi-classical correction, Q, aka 'quantum' correction factor, 
% can be selected from the optional   'type' argument:
%    Q = exp(hw/kT/2)                 'Schofield' or 'Boltzmann'
%    Q = hw/kT./(1-exp(-hw/kT))       'harmonic'  or 'Bader'
%    Q = 2./(1+exp(-hw/kT))           'standard'  or 'Frommhold' (default)
%
% The 'Boltzmann' correction leads to a divergence of the S(q,w) for energies above 
% kT (few 100 meV). The 'harmonic' correction provides a reasonable correction but
% does not fully avoid the divergence at large energies.
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
%   sqw = deBosify(sqw_T, type)
%
% input:
%   sqw_T:  Sqw model (quantum, non classical, including T Bose factor)
%           e.g. 2D model with q as 1st axis (rows, Angs), w as 2nd axis (meV).
%   type: 'Schofield' or 'harmonic' or 'standard' (default)
%
% output:
%   sqw: classical/symmetric Sqw model (iFunc_Sqw2D).
%
% Example: s=sqw_difusion; sb=Bosify(s); sbb=deBosify(sb);
%
% See also: iFunc_Sqw2D, iData_Sqw2D
% (c) E.Farhi, ILL. License: EUPL.

  if nargin < 2, type = 'standard'; end
  
  self = Bosify(self, [ type ' debosify' );
  
  if nargout == 0 && ~isempty(inputname(1)) && isa(self,'iFunc')
    assignin('caller',inputname(1),self);
  end
  
