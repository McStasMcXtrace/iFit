function [Sqw, Iqt, Wq, Tall] = multi_phonons_incoherent(gw, q, T, sigma, m, n)
% multi_phonons_incoherent: compute the multi-phonon contributions in S(q,w) from an initial density of states in the incoherent gaussian approximation
%
% [Sqw, Wq, Tall] = multi_phonons_incoherent(gw, q, T, sigma, m, n)
%
% The partial differential scattering cross section (measured intensity) is:
%
%   (d^2 sigma/dOMEGA dE) = kf/ki S(q,w) 
%
% i.e. the returned S(q,w) is summed with the bound neutron cross section and
% averaged over mode polarisations (e.g. for an isotropic powder).
%
% This implementation is in principle exact for an isotropic monoatomic material,
% e.g. a liquid or powder.
% This methodology is equivalent to the LEAPR module of NJOY ("phonon expansion")
% to compute S(alpha,beta) from a vibrational density of states.
%
% conventions:
% w = Ei-Ef = energy lost by the neutron
%    w > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    w < 0, neutron gains energy, anti-Stokes
%
% Example:
%   s   = sqw_phonons('POSCAR_Al', 'emt');
%   gw  = dos(s);
%   Sqw = multi_phonons_incoherent(gw, [], 300);
%   subplot(Sqw);
%
% input:
%   gw: the vibrational density of states per [meV] [iData]
%   q:  the momentum axis [Angs-1, vector]
%   T:  temperature [K]
%   sigma: neutron cross section [barns]
%   m:  mass [g/mol]
%   n:  number of iterations in the series expansion, e.g. 5
%
% output:
%   Sqw:  neutron weighted S(q,w) terms, to be summed [iData array]
%   Iqt:  I(q,t) terms, to be summed [iData array]
%   Wq:   half Debye-Waller factor. The DW function is exp(-2*Wq)  [iData vs q]
%   Tp:   p-phonon terms [iData array]
%
% Reference:
%   H. Schober, Journal of Neutron Research 17 (2014) 109–357
%     DOI 10.3233/JNR-140016 (see esp. pages 328-331)
%   V.S. Oskotskii, Sov. Phys. Solid State 9 (1967), 420.
%   A. Sjolander, Arkiv for Fysik 14 (1958), 315.

  % check input parameters
  if nargin < 2, q    =[]; end
  if nargin < 3, T    =[]; end
  if nargin < 4, sigma=[]; end
  if nargin < 5, m    =[]; end
  if nargin < 6, n    =[]; end
  
  % search in the vDOS data set in case missing properties are stored there
  % q, T, mass, sigma
  if isempty(q)
    hw = max(abs(gw{1}));
    q  = linspace(0, sqrt(hw), 100);
  end
  if isempty(T) || T<=0
    T = Sqw_getProp(gw);
  end
  if isempty(m) || m<=0
    m = Sqw_getProp(gw, {'Masses','Molar_mass','Mass'});
  end
  if isempty(sigma) || all(sigma<=0)
    sigma = Sqw_getProp(gw, {'sigma_coh','sigma_inc','sigma'});
  end
  % fail when missing information
  if isempty(m)     error([ mfilename ': Unspecified molar mass.' ]); end
  if isempty(T)     error([ mfilename ': Unspecified temperature.' ]); end
  if isempty(sigma) error([ mfilename ': Unspecified scattering cross section.' ]); end
  if isempty(n), n=5; end
  
  % conversion factors
  mn      = 1.674927471E-027; % neutron mass [kg]
  e       = 1.60217662E-019;  % [C]
  HBAR    = 1.05457168e-34;   % Plank/2PI [m2 kg / s]
  kb      = 1.38064852E-023;  % Boltzmann [J/K]
  q2toE   = HBAR*HBAR/2/mn/e*1000*1e20; % = 2.0721 = [Angs^-2] to [meV] 
  meVtoK  = e/1000/kb;        % = 11.6045 = [meV] to [K]
  
  % normalise the input density of states
  % extend the vDOS to the negative side
  titl  = gw.Title;
  gw1   = xlim(gw, [0 inf]);
  gw1   = gw1/trapz(gw1);               % normalize to 1  on [0:inf]
  gw2   = gw1; gw2{1} = -gw2{1};
  gw    = cat(1, gw1, gw2);
  hw    = gw{1};                        % energy axis in [meV]
  
  % we search the smallest non zero energy value and create a 'delta' function
  [dw, index] =  min(abs(hw(:)));
  delta = 0*hw; delta(index) = 1;
  delta = iData(hw, delta);             % make it an iData for easier handling

  % compute the Debye-Waller function W(q) = <u2>/2. Schober p 328 -------------
  kT = T/meVtoK;                        % T [K] -> [meV] = 11.6045
  
  % f(0) is an inverse energy. integrate -inf:inf
  % f(0)=\int(dw g(w)/w [n(w)+1]) Eq (10.58)
  nw = 1./(exp(hw./kT) -1);             % Bose population factor
  % compute the T1 term
  T1 = gw./hw.*(nw+1);                  % Eq. (10.68) p 329
  f0 = trapz(T1);
  disp([ mfilename ': 1/f(0)=' num2str(1/f0) ' [meV]' ]);

  % Debye-Waller coefficient, aka gamma=2W/q^2
  W  = (q2toE/4/m*f0);                  % W = <u2> in Angs^2. integrate 0:inf
  Wq = W.*(q.^2);                       % Eq (10.61), here W(q), not 2W(q). Unit-less
  gamma = 2*W;
  disp([ mfilename ': Debye-Waller coefficient gamma      =' num2str(2*W) ...
    ' [Angs^2]. Debye-Waller function is 2W(q)=gamma*q^2' ]);
  disp([ mfilename ': Debye-Waller factor      exp(-2W(q))=[max mean min]=' ...
    mat2str([max(exp(-2*Wq)) mean(exp(-2*Wq)) min(exp(-2*Wq))],2) ...
    ' on q=[0:' num2str(max(q)) '] [Angs-1]' ]);
  
  Wq = iData(q, Wq);                    % make it an iData for easier handling
  
  % W(q) =<(Q.u)^2> is unit less. <u2>=0.25 in deuterium
  % exp(-2W) is the Debye-Waller factor, unit-less

  % determine the S(q,w) p=0 term in [q,w] space -------------------------------
  Sqw      = (sigma/4/pi)*exp(-2*Wq').*delta';   % eq (10.64) elast.
  Sqw      = Sqw';
  setalias(Sqw, 'Temperature', T, 'Temperature [K]');
  setalias(Sqw, 'Masses',      m, 'Molar masses [g/mol]');
  setalias(Sqw, 'Sigma',   sigma, 'Neutron scattering cross section [barns]');
  Sqw.Title= [ 'Elastic scattering function S(q,w) [p=0] from ' titl ];
  title(Sqw, [ 'S(q,w) [p=0] from ' titl ]);
  Sqw.Label= [ 'S(q,w) [p=0]' ];
  xlabel(Sqw, 'q wavevector [Angs-1]'); ylabel(Sqw, 'hw energy [meV]');
  
  % determine the I(q,t) p=0 term in [q,w] space -------------------------------
  % from T1 we can compute f(t) as the inverse FFT
  % switch T1 energy to Hz so that we get time after the FFT: 1 meV = 241.8 GHz
  if nargout > 1
    T1copy = copyobj(T1);
    setaxis(T1copy, 1, T1copy{1}*241.8e9);  % [meV] -> [Hz]
    ft     = real(fft(T1copy));
    xlabel(ft, 'time [s]');
    % rescale f(t) so that it matches f(0)
    fts    = ft{0};
    ft     = ft/fts(1)*f0;
    % then can compute I(q,t) from Eq (10.63)
    one    = ones(size(ft')); 
    one    = iData(ft{1}, one);             % make it an iData for easier handling
    xlabel(one, 'time [s]');
    Iqt    = sigma/4/pi*exp(-2*Wq').*one;
    Iqt    = Iqt';
    setalias(Iqt, 'Temperature', T, 'Temperature [K]');
    setalias(Iqt, 'Masses',      m, 'Molar masses [g/mol]');
    setalias(Iqt, 'Sigma',   sigma, 'Neutron scattering cross section [barns]');
    Iqt.Title= [ 'Intermediate scattering function I(q,t) [p=0] from ' titl ];
    title(Iqt, [ 'I(q,t) [p=0] from ' titl ]);
    Iqt.Label= [ 'I(q,t) [p=0]' ];
    xlabel(Iqt, 'q wavevector [Angs-1]'); ylabel(Iqt, 'time [s]');
  end
  
  % compute the T[p] terms iteratively -----------------------------------------
  Tall = T1;
  
  % stop when requesting only the one-phonon term
  if isempty(n) || n < 2, return; end
  
  % evaluate iteratively the higher order terms
  Tpm1 = T1;
  fp   = f0.^(1:n);
  for p=2:n
    % compute the T[p] = conv(T[1], T[p-1])
    Tp   = conv(T1, Tpm1);
    % make sure trapz(Tp) = f0^p Eq (10.77)
    Tp   = Tp/trapz(Tp)*fp(p);
    Tpm1 = Tp;
    Tall = [ Tall Tp ];
  end
  
  % compute the S(q,w) expansion [Sjolander] -----------------------------------
  fact  = 1;
  for p=1:n
    fact = fact*p;  % p!
    if numel(Tall) == 1, Tp = Tall; else Tp = Tall(p); end
    
    % we use the normalized T(p) 'tilde' terms Ttilde(p) = T(p)/f(0)^p Eq (10.77)
    Tp      = Tp/fp(p);                 % fp(p) = trapz(Tp) = f0^p Eq (10.78)
    dSqw_q  = exp(-2*Wq).*(2*Wq).^p;    % q-dependence
    dSqw_w  = Tp/trapz(Tp);             % w-dependence Ttilde(p)
    dSqw_c  = sigma/4/pi/fact;          % scaling with constant
    
    % sigma/4/pi*exp(-2*Wq).*(2*Wq).^p / fact * Ttilde_p
    % the (q,w) terms are multiplied orthogonally with iData objects to create a 2D
    
    dSqw    = dSqw_c.*(dSqw_q'.*dSqw_w'); % Eq (10.79) p331 ; [q]*[w]
    dSqw    = dSqw';
    setalias(dSqw, 'Temperature', T, 'Temperature [K]');
    setalias(dSqw, 'Masses',      m, 'Molar masses [g/mol]');
    setalias(dSqw, 'Sigma',   sigma, 'Neutron scattering cross section [barns]');
    dSqw.Title= [ 'Scattering function S(q,w) [p=' num2str(p) ']' ];
    title(dSqw, [ 'S(q,w) [p=' num2str(p) '] from ' titl ]);
    dSqw.Label= [ 'S(q,w) [p=' num2str(p) ']' ];
    xlabel(dSqw, 'q wavevector [Angs-1]'); ylabel(dSqw, 'hw energy [meV]');
    Sqw     = [ Sqw dSqw ];
  end
  
  % compute the I(q,t) expansion [Sjolander] -----------------------------------
  if nargout > 1
    fact= 1;
    for p=1:n
      fact    = fact*p;  % p!
      dIqt_q  = exp(-2*Wq).*(q2toE/2/m.*q.^2).^p; % q-dependence
      dIqt_t  = ft.^p;                            % time dependence
      dIqt_c  = sigma/4/pi/fact;                  % scaling with constant
      dIqt    = dIqt_c.*(dIqt_q'.*dIqt_t');
      setalias(dIqt, 'Temperature', T, 'Temperature [K]');
      setalias(dIqt, 'Masses',      m, 'Molar masses [g/mol]');
      setalias(dIqt, 'Sigma',   sigma, 'Neutron scattering cross section [barns]');
      dIqt = dIqt';
      dIqt.Title= [ 'Intermediate scattering function I(q,t) [p=' num2str(p) ']' ];
      title(dIqt, [ 'I(q,t) [p=' num2str(p) '] from ' titl ]);
      dIqt.Label= [ 'I(q,t) [p=' num2str(p) ']' ];
      xlabel(dIqt, 'q wavevector [Angs-1]'); ylabel(dIqt, 'time [s]');
      Iqt = [ Iqt dIqt ];
    end
  end

% ------------------------------------------------------------------------------
function  T = Sqw_getProp(s, prop)
% Sqw_getProp: search for a property value in a data set
%
% input:
%   s: any iData object, including S(q,w) and DOS ones.
  
  T = [];
  if nargin < 2,    prop = []; end
  if isempty(prop), prop={'Temperature','T'}; end
  if ~iscellstr(prop), prop = { prop }; end
  
  % handle arrays
  if numel(s) > 1
    for index=1:numel(s)
      t = Sqw_getT(s(index), prop);
      if isempty(t), t=nan; end
      T = [ T t ];
    end
    return
  end
  
  for field=prop
    if isfield(s,field{1}), T=get(s,field{1}); end
    
    if isempty(T) || all(T(:)<0) || ~isnumeric(T) || ~isvector(T)
      f = findfield(s,field{1},'exact numeric');
      if ~isempty(f), T = get(s,f{1}); end
    end

    if isempty(T) || all(T(:)<0) || ~isnumeric(T) || ~isvector(T)
      f = findfield(s,field{1},'numeric cache');
      for index=1:numel(f)
        try
          T = get(s,f{1});
          if ~isempty(T) && isnumeric(T) && isvector(T) && all(T(:)>0), break; end
        end
      end
    end

    if ~isempty(T) && isnumeric(T) && isvector(T) && all(T(:)>0), break; end
  end
  
  if isvector(T), T = mean(T(:)); else T=[]; end
  
 
  
% ------------------------------------------------------------------------------
  
% Impulse approx. Schober p 317 
% S(q,w) eq (10.25) Gaussian centred on recoil for Q -> Inf
% W(q)   eq (10.6)  for one single oscillator, isotropic scatterer
%        eq (10.3)  for a xtal

% sin(theta) integrated density of states. Oskotskii/Bredov eq (9.299) p 315
% incoherent approximation for coherent = 20% eq (9.295) p314
%   for an isotropic harmonic oscillator
