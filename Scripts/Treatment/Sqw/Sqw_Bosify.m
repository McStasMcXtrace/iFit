function s = Sqw_Bosify(s, T)
% Sqw_Bosify: apply Bose factor (detailed balance) to a classical data set.
%   The initial data set should obey S(q,w) = S(q,-w), i.e. be 'classical'.
%
% input:
%   s: Sqw data set (classical, symmetric in energy, no T Bose factor)
%        e.g. 2D data set with w as 1st axis (rows), q as 2nd axis.
%   T: when given, Temperature to use for Bose. When not given, the Temperature
%      is searched in the object.
%
% conventions:
% omega = Ei-Ef = energy lost by the neutron
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
% Egelstaff, Eq (9.25) p189
%    S(q,-w) = exp(-hw/kT) S(q,w)
%    S(q,w)  = exp( hw/kT) S(q,-w)
%    S(q,w)  = exp(hw/2kT) S*(q,w) with S*=classical limit
% for omega > 0, S(q,w) > S(q,-w)
%
%  Bose factor: n(w) = 1./(exp(w*11.605/T) -1) ~ exp(-w*11.605/T)
%               w in [meV], T in [K]
%
% Example: s = Sqw_Bosify(s, 300);
%
% See also: Sqw_deBosify, Sqw_symmetrize, Sqw_dynamic_range, Sqw_total_xs

  if nargin == 1, T = []; end

  % handle array of objects
  if numel(s) > 1
    sqw = [];
    for index=1:numel(s)
      sqw = [ sqw feval(mfilename, s(index), T) ];
    end
    s = sqw;
    return
  end

  % check if the data set is Sqw (2D)
  w_present=0;
  q_present=0;
  if isa(s, 'iData') && ndims(s) == 2
    for index=1:2
      lab = lower(label(s,index));
      if any(strfind(lab, 'wavevector')) || any(strfind(lab, 'q')) || any(strfind(lab, 'Angs'))
        q_present=index;
      elseif any(strfind(lab, 'energy')) || any(strfind(lab, 'w')) || any(strfind(lab, 'meV'))
        w_present=index;
      end
    end
  end
  if ~w_present || ~q_present
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not seem to be an isotropic S(|q|,w) 2D object. Ignoring.' ]);
    return
  end
  
  if isempty(T),  T = Sqw_getT(s); end
  if isempty(T) || T == 0
    return
  end

  % test if classical
  if ~isempty(findfield(s, 'classical'))
    if s.classical == 0
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not seem to be classical. It may already contain the Bose factor in which case the detailed balance may be wrong.' ]);
    end
  end

  % check if we need to transpose the S(q,w)
  if w_present==2 && q_present==1
    s = transpose(s);
  end
  
  T2E       = (1/11.6045);           % Kelvin to meV = 1000*K_B/e
  kT        = T*T2E;
  hw_kT     = s{1}/kT;               % hbar omega / kT
  
  % apply sqrt(Bose) factor to get experimental-like
  %n         = exp(hw_kT/2);          % detailed balance (raw)
  n         = hw_kT./(1-exp(-hw_kT));  % Bose factor (true), also satisfies detailed balance = w*(1+n(w))
  n(find(s{1}==0) = 1;
  s         = s .* n;  % apply detailed balance with Bose
  s.Temperature = T;
  s.classical   = 0;
