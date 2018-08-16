function s = Bosify(s0, varargin)
% iData_Sqw2D: Bosify: apply the 'Bose' factor (detailed balance) to a classical data set.
%   The initial data set should obey S*=S(q,w) = S(q,-w), i.e. be 'classical'.
%   The resulting data set is 'quantum/experimental' and satisfies the detailed
%   balance. It contains the temperature effect (population).
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
%   sqw_T = Bosify(sqw)
%   sqw_T = Bosify(sqw, T, type)
%   sqw_T = Bosify(sqw, 'T', T, 'type', type)
%
% input:
%   sqw:  Sqw data set (classical, symmetric in energy, no T Bose factor)
%           e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
%   T:    when given, Temperature to use for Bose. When not given, the Temperature
%           is searched in the object. The temperature is in [K]. 1 meV=11.605 K.
%   type: 'Schofield' or 'harmonic' or 'standard' (default)
%
% output:
%   sqw_T: quantum Sqw data set (non classical, iData_Sqw2D).
%
% Example: s=iData_Sqw2D('SQW_coh_lGe.nc'); sb=Bosify(symmetrize(s), 1235);
%
% See also: iData_Sqw2D/deBosify, iData_Sqw2D/symmetrize, 
%           iData_Sqw2D/dynamic_range, iData_Sqw2D/scattering_cross_section
% (c) E.Farhi, ILL. License: EUPL.

% classical   T     getT    can convert
% --------------------------------------------------------------------
% yes         []    any     yes, T=getT would always be >0
% yes         T             yes, use any T <0 or >0
% no(quantum) []    any     must use T=-getT to remove
% no(quantum) T     any     must use -abs(T) to remove
% no(quantum) T>0           NO
  
  
  p = varargin2struct({'T' 'type'}, varargin, true);
  if isempty(p.type), p.type='standard'; end
  
  s=[];
  
  % handle array of objects
  if numel(s0) > 1
    for index=1:numel(s)
      s = [ s feval(mfilename, s0(index), p.t, p.type) ];
    end
    return
  end
  
  if isempty(s0), return; end
  s = copyobj(s0);
  if strfind(lower(p.type),'debosify')
       do_bosify=0;
       p.type = strtrim(strrep(lower(p.type), 'debosify','')); % remove debosify occurence
  else do_bosify=1; end

  % define the fallback Temperature to use
  if isempty(p.t) || p.t == 0
    T0 = Sqw_getT(s); 
    if (isempty(T0) || T0 == 0)
      T0 = 293;
    end
    disp([ mfilename ': INFO: Using Temperature=' num2str(T0) ' [K] for data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
    p.t=T0;
  else T0=p.t; end  
  
  % when not given use the p.t from the data set.
  if isempty(p.t) || p.t == 0
    error([ mfilename ': ERROR: Undefined Temperature in data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
  end

  % test if classical/quantum and bosify/debosify agree
  if isfield(s0,'classical') || ~isempty(findfield(s0, 'classical'))
    % Bosify   must be applied on classical
    % deBosify must be applied on quantum
    classical = get(s0,'classical');
    if (~isempty(classical) && classical(1) == 0 && do_bosify)
      disp([ mfilename ': WARNING: Not "classical/symmetric": The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not seem to be classical (classical=0).' ]);
      disp([ mfilename ':   It may ALREADY contain the Bose factor in which case the detailed balance will be wrong.' ]);
    elseif (~isempty(classical) && classical(1) == 1 && ~do_bosify)
      disp([ 'de' mfilename ': WARNING: Not "quantum": The data set ' s.Tag ' ' s.Title ' from ' s.Source ' seems to be classical/symmetric (classical=1).' ]);
      disp([ 'de' mfilename ':   The Bose factor may NOT NEED to be removed in which case the detailed balance will be wrong.' ]);
    end
  end
  
  % handle different temperatures
  if numel(p.t) > 1
    sqw = [];
    for index=1:numel(p.t)
      sqw = [ sqw feval(mfilename, s, p.t(index), p.type, 'checked') ];
    end
    s = sqw;
    return
  end
  
  T2E       = (1/11.6045);           % Kelvin to meV = 1000*K_B/e
  kT        = p.t*T2E;
  hw_kT     = getaxis(s,1)./kT;               % hbar omega / kT
  
  % apply sqrt(Bose) factor to get experimental-like
  % semi-classical corrections, aka quantum correction factor
  
  if strcmpi(p.type, 'Schofield') || strcmpi(p.type, 'exp') || strcmpi(p.type, 'Boltzmann') 
    % P. Schofield. Phys. Rev. Lett., 4, 239 (1960).
    Q  = exp(hw_kT/2);          
  elseif strcmpi(p.type, 'Harmonic') || strcmpi(p.type, 'Bader')
    % J. S. Bader and B. J. Berne. J. Chem. Phys., 100, 8359 (1994).
    % p.t. D. Hone and G. A. Voth. J. Chem. Phys., 121, 6412 (2004).
    Q  = hw_kT./(1-exp(-hw_kT));  % w*(1+n(w))
  elseif strcmpi(p.type, 'Standard') || strcmpi(p.type, 'Frommhold') || strcmpi(p.type, 'default')
    % L. Frommhold. Collision-induced absorption in gases, 1 st ed., Cambridge
    %   Monographs on Atomic, Molecular, and Chemical Physics, Vol. 2,
    %   Cambridge Univ. Press: London (1993).
    Q = 2./(1+exp(-hw_kT));
  else
    error([ mfilename ': Unknown semi-classical correction ' p.type ]);
  end
  Q(find(hw_kT==0)) = 1;
  if ~do_bosify, Q = 1./Q; end
  s         = s .* Q;  % apply detailed balance with the selected correction
  
  setalias(s, 'Temperature', abs(p.t));
  setalias(s, 'QuantumCorrection',p.type,[ 'Quantum correction applied in ' mfilename ]);
  if do_bosify % Bosify
    setalias(s, 'classical', 0, 'This is a experimental/quantum S(q,w) with Bose factor');
    s = commandhistory(s, sprintf('Bosify(%s,%g,''%s'')', s0.Tag,p.t,p.type));
  elseif ~do_bosify % deBosify
    setalias(s, 'classical', 1, 'This is a classical/symmetric S(q,w) without Bose factor');
    s = commandhistory(s, sprintf('deBosify(%s,%g,''%s'')', s0.Tag, p.t,p.type));
  end
  
  s = iData_Sqw2D(s); % final Sqw2D object
