function s = Bosify(s0, T, type)
% Bosify: apply the 'Bose' factor (detailed balance) to a classical data set.
%
% iData_Sqw2D: bosify: apply the Bose factor, which adds a temperature effect.
%   The initial data set should obey S*=S(q,w) = S(q,-w), i.e. be 'classical'.
%   The S(q,w) is a dynamic structure factor aka scattering function.
%
% sb=Bosify(s, T, type) 
%
% input:
%   s: Sqw data set (classical, symmetric in energy, no T Bose factor)
%        e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
%   T: when given, Temperature to use for Bose. When not given, the Temperature
%      is searched in the object. The temperature is in [K]. 1 meV=11.605 K.
%   type: 'Schofield' or 'harmonic' or 'standard' (default)
%
% output:
%   sb: quantum Sqw data set (non classical, iData_Sqw2D).
%
% conventions:
% omega = Ei-Ef = energy lost by the neutron, given in [meV]
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
% Egelstaff, Eq (9.25) p189,  where Q(w) is defined below:
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
% Example: s=iData_Sqw2D('SQW_coh_lGe.nc'); sb=Bosify(s, 300);
%
% See also: deBosify, symmetrize, Sqw_dynamic_range, Sqw_scatt_xs
% (c) E.Farhi, ILL. License: EUPL.

% classical   T     getT    can convert
% --------------------------------------------------------------------
% yes         []    any     yes, T=getT would always be >0
% yes         T             yes, use any T <0 or >0
% no(quantum) []    any     must use T=-getT to remove
% no(quantum) T     any     must use -abs(T) to remove
% no(quantum) T>0           NO

  if nargin < 2, T = []; end
  if nargin < 3, type=''; end
  if isempty(type), type='standard'; end
  s=[];
  
  % handle array of objects
  if numel(s0) > 1
    for index=1:numel(s)
      s = [ s feval(mfilename, s0(index), T, type) ];
    end
    return
  end
  
  if isempty(s0), return; end
  s = copyobj(iData(s0)); % back to iData

% define the fallback Temperature to use
  if isempty(T),    T0 = Sqw_getT(s); else T0=T; end  
  if isempty(T) && (isempty(T0) || T0 == 0)
    T0 = 300;
  end
  
  % adapt the Temperature to the type of classical/quantum
  if ~isempty(T0) && ...
    (isfield(s,'classical') || ~isempty(findfield(s, 'classical')))
    if     s.classical == 1,     T0= abs(T0);
    elseif s.classical == 0,     T0=-abs(T0); end
  end
  
  % when not given use the T from the data set.
  if isempty(T)
    T=T0; 
    if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
      if s.classical, t='classical'; else t='quantum'; end
      if T0<0, d='de'; else d=''; end
      disp([ d mfilename ': WARNING: Using Temperature=' num2str(abs(T0)) ' [K] for ' t ' data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
    end
  end
  if isempty(T) || T == 0
    error([ mfilename ': WARNING: Undefined Temperature in data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
  end

  % test if classical and T agree: 
  if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
    if (s.classical == 0 && T > 0)
      disp([ mfilename ': WARNING: Not "classical/symmetric": The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not seem to be classical (classical=0).' ]);
      disp([ mfilename ':   It may ALREADY contain the Bose factor in which case the detailed balance will be wrong.' ]);
    elseif (s.classical == 1 && T < 0)
      disp([ 'de' mfilename ': WARNING: Not "quantum": The data set ' s.Tag ' ' s.Title ' from ' s.Source ' seems to be classical/symmetric (classical=1).' ]);
      disp([ 'de' mfilename ':   The Bose factor may NOT NEED to be removed in which case the detailed balance will be wrong.' ]);
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
  
  setalias(s, 'Temperature', abs(T));
  setalias(s, 'QuantumCorrection',type,[ 'Quantum correction applied in ' mfilename ]);
  if T> 0 % Bosify
    setalias(s, 'classical', 0, 'This is a quantum S(q,w) with Bose factor');
    s = commandhistory(s, sprintf('Bosify(%s,%g,''%s'')', s0.Tag,T,type));
  elseif T< 0 % deBosify
    setalias(s, 'classical', 1, 'This is a classical/symmetric S(q,w) without Bose factor');
    s = commandhistory(s, sprintf('deBosify(%s,%g,''%s'')', s0.Tag, -T,type));
  end
  
  s = iData_Sqw2D(s); % final Sqw2D object
