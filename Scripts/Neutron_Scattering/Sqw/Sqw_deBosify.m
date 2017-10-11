function s = Sqw_deBosify(s, T, type)
% Sqw_deBosify: remove Bose factor (detailed balance) from an 'experimental' data set.
%  In principle the resulting data set is 'classical' that is S(q,w) = S(q,-w)
%
%  The S(q,w) is a dynamic structure factor aka scattering function.
%
% input:
%   s: Sqw data set (non classical, including T Bose factor e.g from experiment)
%        e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
%   T: when given, Temperature to use for Bose. When not given, the Temperature
%      is searched in the object.
%   type: 'Schofield' or 'harmonic' or 'standard' (default)
%
% conventions:
% omega = Ei-Ef = energy lost by the neutron, given in [meV]
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
% Egelstaff, Eq (9.25) p189
%    S(q,-w) = exp(-hw/kT) S(q,w)
%    S(q,w)  = exp( hw/kT) S(q,-w)
%    S(q,w)  = Q(w) S*(q,w) with S*=classical limit
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
% Example: s = Sqw_deBosify(s, 300);
%
% See also: Sqw_Bosify, Sqw_symmetrize, Sqw_dynamic_range, Sqw_scatt_xs
% (c) E.Farhi, ILL. License: EUPL.

  if nargin == 0, return; end
  if ~isa(s, 'iData'), s=iData(s); end
  
  if nargin < 2, T = []; end
  if nargin < 3, type=''; end
  if isempty(type), type='standard'; end

  % handle array of objects
  if numel(s) > 1
    sqw = [];
    for index=1:numel(s)
      sqw = [ sqw feval(mfilename, s(index), T, type) ];
    end
    s = sqw;
    return
  end

  s = Sqw_check(s);
  if isempty(s), return; end
  
  if isempty(T) || T == 0
    return
  end
  
  % test if classical
  if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
    if s.classical == 1
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' seems to already be classical. The detailed balance removal may be wrong.' ]);
    end
  end
  
  % get symmetric from experimental data
  s = Sqw_Bosify(s, -T, type);
  setalias(s, 'Temperature', T);
  setalias(s, 'QuantumCorrection',type,[ 'Quantum correction applied in ' mfilename ]);
  setalias(s, 'classical', 1);
  
