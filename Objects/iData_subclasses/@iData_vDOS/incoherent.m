function [Sqw, Iqt, Wq, Tall] = incoherent(gw, varargin)
% iData_vDOS: incoherent: compute the incoherent gaussian approximation scattering law S(q,w) from an initial density of states vDOS
%
%   [Sqw, Iqt, Wq, Tall] = incoherent(gw, q, T, m, n, DW)
%
% compute: Density of States -> Incoherent approximation S(q,w)
%
% The input argument 'gw' should be a vibrational density of states (vDOS) as
%   obtained from an experiment (e.g. Bedov/Oskotskii estimate), molecular 
%   dynamics, or lattice dynamics. In the so-called incoherent approximation,
%   the vDOS obtained from incoherent and coherent scattering laws are equal.
%
% The result is the dynamic structure factor (scattering law) for neutrons, in
%   the incoherent gaussian approximation. The corresponding intermediate 
%   scattering function is also returned. These should be e.g. multiplied by the
%   neutron scattering bound cross section 'sigma_inc' [barns]. This calculation
%   includes the Debye-Waller factor.
%
% This implementation is in principle exact for an isotropic monoatomic material,
%   e.g. a liquid, powder, or cubic crystal. 
% This methodology is equivalent to the LEAPR module of NJOY ("phonon expansion")
%   to compute S(alpha,beta) from a vibrational density of states.
%
% For a poly-atomic material with a set of non-equivalent atoms with relative 
%   concentration Ci and mass Mi and bound scattering cross section sigma_i, 
%   one should use:
%
%   m     = [sum_i Ci sigma_i]/[sum_i Ci sigma_i/Mi]      weighted mass
%   sigma = sum_i Ci sigma_i                              weighted cross section
%
% conventions:
% w = Ei-Ef = energy lost by the neutron
%    w > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    w < 0, neutron gains energy, anti-Stokes
%
% Theory:
%   T1(w) = g(w)./hw.*(nw+1)
%   f0    = \int T1(w) dw
%   fp    = f0^p
%   W(q)  = h^2*q^2/2/m*f(0)
%   Tp    = conv( T1, Tp-1 )
%
% The first term returned is the Elastic Incoherent [p=0]
%    S(q,w)[p=0] = (1/4/pi)*exp(-2*W(q)).*delta(w)
% the other terms are obtained by iterative auto-convolution by the vDOS
%    S(q,w)[p]   = 1/4pi/!p exp(-2*W(q)).*(2*W(q)).^p*Tp(w)
%
% Reference:
%   H. Schober, Journal of Neutron Research 17 (2014) 109–357
%     DOI 10.3233/JNR-140016 (see esp. pages 328-331)
%   V.S. Oskotskii, Sov. Phys. Solid State 9 (1967), 420.
%   A. Sjolander, Arkiv for Fysik 14 (1958), 315.
%
% syntax:
%   [Sqw, Iqt, Wq, Tall] = incoherent(gw)
%   [Sqw, Iqt, Wq, Tall] = incoherent(gw, q, T, m, n, DW)
%   [Sqw, Iqt, Wq, Tall] = incoherent(gw, 'q', q, 'T', T, 'm', m, 'n', n, 'DW', dw)
%
% Missing arguments (or given as [] empty), are searched within the initial density 
%   of states object. Input arguments can be given in order, or with name-value 
%   pairs, or as a structure with named fields.
%
% input:
%   gw: the vibrational density of states per [meV] [iData]
%   q:  the momentum axis [Angs-1, vector]
%   T:  temperature [K]
%   m:  mass of the scattering unit [g/mol]
%   n:  number of iterations in the series expansion, e.g. 5
%   DW: Debye-Waller coefficient gamma=<u^2> [Angs^2] e.g. 0.005
%       The Debye-Waller function is      2W(q)=gamma*q^2
%       The Debye-Waller factor   is exp(-2W(q))
%
% output:
%   Sqw:  neutron weighted S(q,w) incoherent terms, to be summed [iData_Sqw2D array]
%   Iqt:  I(q,t) incoherent terms, to be summed [iData array]
%   Wq:   half Debye-Waller factor. The DW function is exp(-2*Wq)  [iData vs q]
%   Tp:   p-phonon terms [iData array]
%
% Example:
%   s   = sqw_phonons('POSCAR_Al', 'emt');
%   gw  = dos(s);
%   Sqw = incoherent(gw, [], 300);
%   subplot(Sqw);
%
% See also: iData_Sqw2D/multi_phonons
% (c) E.Farhi, ILL. License: EUPL.
  
  Sqw=[]; Iqt=[]; Wq=[]; Tall=[];
  if isempty(gw), return; end
  pars = varargin2struct({'q' 'T' 'm' 'n' 'DW'}, varargin, true);
  
  if isfield(pars, 'gamma') && ~isempty(pars.gamma), pars.dw = pars.gamma; end
  if isfield(pars, 'u2')    && ~isempty(pars.u2),    pars.dw = pars.u2; end
  
  % search in the vDOS data set in case missing properties are stored there
  % q, T, mass
  if isempty(pars.q)
    hw = max(abs(getaxis(gw,1)));
    pars.q  = linspace(0, sqrt(hw), 100);
  end
  if isempty(pars.t) || pars.t<=0
    pars.t = Sqw_getT(gw);
  end
  if isempty(pars.m) || pars.m<=0
    pars.m = sum(Sqw_getT(gw, {'Molar_mass','Masses','Mass','Weight'}));
  end
 
  % fail when missing information
  if isempty(pars.m) || pars.m<=0    error([ mfilename ': Unspecified molar mass (m). Use incoherent(g, q, T, m)' ]); end
  if isempty(pars.t) || pars.t<=0    error([ mfilename ': Unspecified temperature (T). Use incoherent(gw, q, T).' ]); end
  if isempty(pars.n),   pars.n=5; end
  
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
  gw2   = gw1; setaxis(gw2,1, -getaxis(gw2,1));
  gw    = cat(1, gw1, gw2);
  hw    = getaxis(gw,1);                % energy axis in [meV]
  
  % we search the smallest non zero energy value and create a 'delta' function
  [dw, index] =  min(abs(hw(:)));
  delta = 0*hw; delta(index) = 1/dw;    % delta has a norm = 1
  delta = iData(hw, delta);             % make it an iData for easier handling

  % compute the Debye-Waller function W(q) = <u2>/2. Schober p 328 -------------
  kT = pars.t/meVtoK;                        % T [K] -> [meV] = 11.6045
  
  % f(0) is an inverse energy. integrate -inf:inf
  % f(0)=\int(dw g(w)/w [pars.n(w)+1]) Eq (10.58)
  nw = 1./(exp(hw./kT) -1);             % Bose population factor
  % compute the T1 term
  T1 = gw./hw.*(nw+1);                  % Eq. (10.68) p 329
  S.type='()';
  S.subs={ ~isfinite(T1) };
  T1 = subsasgn(T1, S, 0);
  f0 = trapz(T1);
  if ~isempty(pars.dw) && pars.dw > 0
    % we use a specified DW 'gamma' value
    T1 = T1/f0;
    f0 = pars.dw/(q2toE/4/pars.m)/2;  % new DW factor sets f0=trapz(T1)
    T1 = T1*f0;
  end
  disp([ mfilename ': 1/f(0)=' num2str(1/f0) ' [meV]' ]);
  % compute the single oscillator frequency
  % eq (10.62): W    = h2/2/m f(0) = (q2toE/4/pars.m)*f0
  %             W(q) = W q^2
  % 1/f0 = w0/coth(hw0 / 2kT) = w0/[2 n(w0) +1 ]
  % we search for index so that it is satisfied
  [~, index] = min( abs(1/f0 -hw./(2*nw+1)) );
  w0         = abs(hw(index));
  disp([ mfilename ':    w0 =' num2str(w0) ' [meV] single harmonic oscillator frequency' ]);

  % Debye-Waller coefficient, aka gamma=2W/q^2
  W  = (q2toE/4/pars.m*f0);                  % W = <u2> in Angs^2. integrate 0:inf
  Wq = W.*(pars.q.^2);                       % Eq (10.61), here W(q), not 2W(q). Unit-less
  gamma = 2*W;
  disp([ mfilename ': Debye-Waller coefficient gamma=<u^2>=' num2str(2*W) ...
    ' [Angs^2]. Debye-Waller function is 2W(q)=gamma*q^2' ]);
  % gamma=0.00677 Angs^2 for Vanadium (Skold PRA 1972)
  disp([ mfilename ': Debye-Waller factor      exp(-2W(q))=[max mean min]=' ...
    mat2str([max(exp(-2*Wq)) mean(exp(-2*Wq)) min(exp(-2*Wq))],2) ...
    ' on q=[0:' num2str(max(pars.q)) '] [Angs-1]' ]);
  
  Wq = iData(pars.q, Wq);                    % make it an iData for easier handling
  
  % W(q) =<(Q.u)^2> is unit less. <u2>=0.25 in deuterium
  % exp(-2W) is the Debye-Waller factor, unit-less
  
  % in the following, we assume sigma_inc = 1. Result should be scaled.

  % determine the S(q,w) p=0 term in [q,w] space -------------------------------
  Sqw      = (1/4/pi)*exp(-2*Wq').*delta';   % eq (10.64) elast.
  Sqw      = Sqw';
  setalias(Sqw, 'Temperature', pars.t, 'Temperature [K]');
  setalias(Sqw, 'Weight',      pars.m, 'Molar masses [g/mol]');
  setalias(Sqw, 'Molar_mass',  pars.m, 'Molar masses [g/mol]');
  setalias(Sqw, 'DebyeWallerCoeficient', gamma, ...
      'Debye-Waller coefficient gamma=<u^2> [Angs^2]. Debye-Waller function is 2W(q)=gamma*q^2. Debye-Waller factor is exp(-gamma.q^2)');
  Sqw.Title= [ 'Elastic scattering function S(q,w) [p=0] from ' titl ];
  title(Sqw, [ 'S(q,w) [p=0] from ' titl ]);
  Sqw.Label= [ 'S(q,w) [p=0]' ];
  xlabel(Sqw, 'wavevector [Angs-1]'); ylabel(Sqw, 'energy [meV]');
  Sqw      = iData_Sqw2D(Sqw);
  
  % transfer UserData stuff
  for f={'properties','calc','configuration','options','FORCES','dir','maxFreq'}
    if isfield(gw.UserData, f{1})
      Sqw.UserData.(f{1}) = gw.UserData.(f{1});
    end
  end
  
  % determine the I(q,t) p=0 term in [q,w] space -------------------------------
  % from T1 we can compute f(t) as the inverse FFT
  % switch T1 energy to Hz so that we get time after the FFT: 1 meV = 241.8 GHz
  if nargout > 1
    T1copy = copyobj(T1);
    setaxis(T1copy, 1, T1copy{1}*241.8e9);  % [meV] -> [Hz]
    ft     = real(fft(T1copy));
    xlabel(ft, 'time [s]');
    % rescale f(t) so that it matches f(0)
    fts    = getaxis(ft,0);
    ft     = ft/fts(1)*f0;
    % then can compute I(q,t) from Eq (10.63)
    one    = ones(size(ft')); 
    one    = iData(getaxis(ft,1), one);             % make it an iData for easier handling
    xlabel(one, 'time [s]');
    Iqt    = 1/4/pi*exp(-2*Wq').*one;
    Iqt    = Iqt';
    setalias(Iqt, 'Temperature', pars.t, 'Temperature [K]');
    setalias(Iqt, 'Weight',      pars.m, 'Molar masses [g/mol]');
    setalias(Iqt, 'DebyeWallerCoeficient', gamma, ...
      'Debye-Waller coefficient gamma=<u^2> [Angs^2]. Debye-Waller function is 2W(q)=gamma*q^2. Debye-Waller factor is exp(-gamma.q^2)');
    Iqt.Title= [ 'Intermediate scattering function I(q,t) [p=0] from ' titl ];
    title(Iqt, [ 'I(q,t) [p=0] from ' titl ]);
    Iqt.Label= [ 'I(q,t) [p=0]' ];
    xlabel(Iqt, 'q wavevector [Angs-1]'); ylabel(Iqt, 'time [s]');
  end
  
  % compute the T[p] terms iteratively -----------------------------------------
  Tall = T1;
  
  % stop when requesting only the one-phonon term
  if isempty(pars.n) || pars.n < 2, return; end
  
  % evaluate iteratively the higher order terms
  Tpm1 = T1;
  fp   = f0.^(1:pars.n);
  for p=2:pars.n
    % compute the T[p] = conv(T[1], T[p-1])
    Tp   = conv(T1, Tpm1);
    % make sure trapz(Tp) = f0^p Eq (10.77)
    Tp   = Tp/trapz(Tp)*fp(p);
    Tpm1 = Tp;
    Tall = [ Tall Tp ];
  end
  
  % compute the S(q,w) expansion [Sjolander] -----------------------------------
  fact  = 1;
  for p=1:pars.n
    fact = fact*p;  % p!
    if numel(Tall) == 1, Tp = Tall; else Tp = Tall(p); end
    
    % we use the normalized T(p) 'tilde' terms Ttilde(p) = T(p)/f(0)^p Eq (10.77)
    Tp      = Tp/fp(p);                 % fp(p) = trapz(Tp) = f0^p Eq (10.78)
    dSqw_q  = exp(-2*Wq).*(2*Wq).^p;    % q-dependence
    dSqw_w  = Tp/trapz(Tp);             % w-dependence Ttilde(p)
    dSqw_c  = 1/4/pi/fact;     % scaling with constant
    
    % sigma/4/pi*exp(-2*Wq).*(2*Wq).^p / fact * Ttilde_p
    % the (q,w) terms are multiplied orthogonally with iData objects to create a 2D
    
    dSqw    = dSqw_c.*(dSqw_q'.*dSqw_w'); % Eq (10.79) p331 ; [q]*[w]
    dSqw    = dSqw';
    setalias(dSqw, 'Temperature', pars.t, 'Temperature [K]');
    setalias(dSqw, 'Weight',      pars.m, 'Molar masses [g/mol]');
    setalias(dSqw, 'DebyeWallerCoeficient', gamma, ...
      'Debye-Waller coefficient gamma=<u^2> [Angs^2]. Debye-Waller function is 2W(q)=gamma*q^2. Debye-Waller factor is exp(-gamma.q^2)');
    dSqw.Title= [ 'Scattering function S(q,w) [p=' num2str(p) ']' ];
    title(dSqw, [ 'S(q,w) [p=' num2str(p) '] from ' titl ]);
    dSqw.Label= [ 'S(q,w) [p=' num2str(p) ']' ];
    xlabel(dSqw, 'q wavevector [Angs-1]'); ylabel(dSqw, 'hw energy [meV]');
    dSqw    = iData_Sqw2D(dSqw);
    Sqw     = [ Sqw dSqw ];
  end
  
  % compute the I(q,t) expansion [Sjolander] -----------------------------------
  if nargout > 1
    fact= 1;
    for p=1:pars.n
      fact    = fact*p;  % p!
      dIqt_q  = exp(-2*Wq).*(q2toE/2/pars.m.*pars.q.^2).^p; % q-dependence
      dIqt_t  = ft.^p;                            % time dependence
      dIqt_c  = 1/4/pi/fact;                  % scaling with constant
      dIqt    = dIqt_c.*(dIqt_q'.*dIqt_t');
      setalias(dIqt, 'Temperature', pars.t, 'Temperature [K]');
      setalias(dIqt, 'Weight',      pars.m, 'Molar masses [g/mol]');
      setalias(dIqt, 'DebyeWallerCoeficient', gamma, ...
      'Debye-Waller coefficient gamma=<u^2> [Angs^2]. Debye-Waller function is 2W(q)=gamma*q^2. Debye-Waller factor is exp(-gamma.q^2)');
      dIqt = dIqt';
      dIqt.Title= [ 'Intermediate scattering function I(q,t) [p=' num2str(p) ']' ];
      title(dIqt, [ 'I(q,t) [p=' num2str(p) '] from ' titl ]);
      dIqt.Label= [ 'I(q,t) [p=' num2str(p) ']' ];
      xlabel(dIqt, 'q wavevector [Angs-1]'); ylabel(dIqt, 'time [s]');
      Iqt = [ Iqt dIqt ];
    end
  end
  
  if nargout == 0
    fig=figure; 
    G1 =plus(Sqw); title(G1, ['Sinc [single+multi] from ' title(gw) ]); G1.Label=title(G1); G1.Title=title(G1);
    h  =plot(log10(G1)); 
    set(fig, 'NextPlot','new');
  end
  
% ------------------------------------------------------------------------------
  
% Impulse approx. Schober p 317 
% S(q,w) eq (10.25) Gaussian centred on recoil for Q -> Inf
% W(q)   eq (10.6)  for one single oscillator, isotropic scatterer
%        eq (10.3)  for a xtal

% sin(theta) integrated density of states. Oskotskii/Bredov eq (9.299) p 315
% incoherent approximation for coherent = 20% eq (9.295) p314
%   for an isotropic harmonic oscillator
