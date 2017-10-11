function s = Sqw_Bosify(s, T, type, options)
% Sqw_Bosify: apply the 'Bose' factor (detailed balance) to a classical data set.
%   The initial data set should obey S(q,w) = S(q,-w), i.e. be 'classical'.
%
%  The S(q,w) is a dynamic structure factor aka scattering function.
%
%  Sqw_Bosify(Sqw, 0) makes a check of the input data set and returns.
%
% input:
%   s: Sqw data set (classical, symmetric in energy, no T Bose factor)
%        e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
%   T: when given, Temperature to use for Bose. When not given, the Temperature
%      is searched in the object. The temperature is in [K]. 1 meV=11.605 K.
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
% Example: s = Sqw_Bosify(s, 300);
%
% See also: Sqw_deBosify, Sqw_symmetrize, Sqw_dynamic_range, Sqw_scatt_xs
% (c) E.Farhi, ILL. License: EUPL.

  if nargin == 0, return; end
  if ~isa(s, 'iData'), s=iData(s); end
  if nargin < 2, T = []; end
  if nargin < 3, type=''; end
  if nargin < 4, options=''; end

  % handle array of objects
  if numel(s) > 1
    sqw = [];
    for index=1:numel(s)
      sqw = [ sqw feval(mfilename, s(index), T, type) ];
    end
    s(index)=iData; % free memory
    s = sqw;
    return
  end

  if ~strcmp(options, 'checked')
    s = Sqw_check(s); % in private
  end

  if isempty(s), return; end
  
  if isempty(T),  T = Sqw_getT(s); end
  if isempty(type), type='standard'; end
  
  if isempty(T) || T == 0
    return
  end

  % test if classical
  if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
    if s.classical == 0
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not seem to be classical. It may already contain the Bose factor in which case the detailed balance may be wrong.' ]);
    end
  end
  
  % handle different temperatures
  if numel(T) > 1
    sqw = [];
    for index=1:numel(T)
      sqw = [ sqw feval(mfilename, s, T(index), type, 'checked') ];
    end
    s = sqw;
    return
  end
  
  T2E       = (1/11.6045);           % Kelvin to meV = 1000*K_B/e
  kT        = T*T2E;
  hw_kT     = s{1}./kT;               % hbar omega / kT
  
  % apply sqrt(Bose) factor to get experimental-like
  % semi-classical corrections, aka quantum correction factor
  
  if strcmpi(type, 'Schofield') || strcmpi(type, 'exp') || strcmpi(type, 'Boltzmann') 
    % P. Schofield. Phys. Rev. Lett., 4, 239 (1960).
    Q  = exp(hw_kT/2);          
  elseif strcmpi(type, 'harmonic') || strcmpi(type, 'Bader')
    % J. S. Bader and B. J. Berne. J. Chem. Phys., 100, 8359 (1994).
    % T. D. Hone and G. A. Voth. J. Chem. Phys., 121, 6412 (2004).
    Q  = hw_kT./(1-exp(-hw_kT));  % w*(1+n(w))
  elseif strcmpi(type, 'standard') || strcmpi(type, 'Frommhold')
    % L. Frommhold. Collision-induced absorption in gases, 1 st ed., Cambridge
    %   Monographs on Atomic, Molecular, and Chemical Physics, Vol. 2,
    %   Cambridge Univ. Press: London (1993).
    Q = 2./(1+exp(-hw_kT));
  else
    error([ mfilename 'Unknown semi-classical correction ' type ]);
  end
  Q(find(s{1}==0)) = 1;
  s         = s .* Q;  % apply detailed balance with the selected correction
  
  setalias(s, 'Temperature', T);
  setalias(s, 'QuantumCorrection',type,[ 'Quantum correction applied in ' mfilename ]);
  setalias(s, 'classical', 0);
