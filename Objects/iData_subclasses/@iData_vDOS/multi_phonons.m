function [Gw, Tsym] = multi_phonons(gw, varargin)
% iData_vDOS: multi_phonons: compute the integrated multi-phonon generalised density of states (gDOS) from an initial vibrational density of states (vDOS)
%
%   Gw = multi_phonons(gw, Ki, T, m, n)
%
% compute: Density of States -> multi-phonon gDOS terms
%
% The input argument 'gw' should be a vibrational density of states (vDOS) as
%   obtained from an experiment (e.g. Bedov/Oskotskii estimate), molecular 
%   dynamics, or lattice dynamics. In the so-called incoherent approximation,
%   the vDOS obtained from incoherent and coherent scattering laws are equal.
% This method is equivalent to the 'incoherent' one, but provides the
%   density-of-states gDOS instead of a scattering law Sinc(q,w).
%
% The result is the 'neutron weighted' (generalised) density of states (gDOS),
%   which should be compared with:
%
%   Gw -> \int q. I(q,w) /Ki/Kf dq ~ \int I(theta, w) sin(theta) d(theta)
%
%   where I(q,w) and I(theta,w) are the measured scattering intensity as a function 
%   of the momentum q or scattering angle theta.
%   In practice, the measured I(q,w) contains both the single phonon and multi-phonon
%   terms (and potentially multiple scattering).
%   These should be e.g. multiplied by the neutron scattering bound cross section
%   'sigma' [barns]. This calculation includes the Debye-Waller factor.
%
% The 1st element of the result 'Gw' is the generalised density of states, as computed
%   in the Oskotskii formalism. The following terms are the multi-phonon contributions.
%   The generalised density of states (gDOS) is the sum of the terms in this expansion.
%
% This implementation is in principle exact for an isotropic monoatomic material,
%   e.g. a liquid, powder, or cubic crystal. This methodology is equivalent to 
%   an iteration of the MUPHOCOR code.
%
% For a poly-atomic material with a set of non-equivalent atoms with relative 
%   concentration Ci, mass Mi and bound scattering cross section sigma_i, 
%   one should use:
%
%   sigma = sum_i Ci sigma_i                              weighted cross section
%   m     = [sum_i Ci sigma_i]/[sum_i Ci sigma_i/Mi]      weighted mass
%
% Reference:
%   H. Schober, Journal of Neutron Research 17 (2014) 109–357
%     DOI 10.3233/JNR-140016 (see esp. pages 328-331)
%   V.S. Oskotskii, Sov. Phys. Solid State 9 (1967), 420.
%   A. Sjolander, Arkiv for Fysik 14 (1958), 315.
%   W. Reichardt, MUPHOCOR Karlsruhe Report 13.03.01p06L (1984)
%
% syntax:
%   Gw = multi_phonons(gw, Ki, T, m, n)
%   Gw = multi_phonons(gw, 'Ki',Ki, 'T',T, 'm',m, 'n',n)
%   Gw = multi_phonons(gw, 'lambda', lambda)
%
% Missing arguments (or given as [] empty), are searched within the initial density 
%   of states object. Input arguments can be given in order, or with name-value 
%   pairs, or as a structure with named fields.
%
% input:
%   gw:   the vibrational density of states per [meV] [iData]
%   Ki:   incident wavevector [Angs-1]. Ki=2*pi/lambda=0.695*sqrt(Ei)
%   T:    temperature [K]
%   m:    mass [g/mol]
%   phi:  detector angles, min and max [deg]
%   n:    number of iterations in the series expansion, e.g. 5
%   DW: Debye-Waller coefficient gamma=<u^2> [Angs^2] e.g. 0.005
%       The Debye-Waller function is      2W(q)=gamma*q^2
%       The Debye-Waller factor   is exp(-2W(q))
%   'lambda','Ei': additional named arguments to specify the incident energy
%
% output:
%   Gw:   neutron weighted gDOS terms, to be summed [iData_vDOS array]
%   Tp:   p-phonon terms [iData array]
%
% Example: s=sqw_cubic_monoatomic; gw=dos(s); 
%   mp=multi_phonons(gw); gdos=plus(mp); plot( [gw mp ]);
% 
%
% See also: iData_Sqw2D/dos, iData_vDOS/incoherent
% (c) E.Farhi, ILL. License: EUPL.

  Gw = []; Tsym=[];
  if isempty(gw), return; end
  
  % check input parameters
  pars = varargin2struct({'Ki' 'T' 'm' 'phi' 'n' 'DW' 'lambda' 'Ei'}, varargin, true);
  
  if isfield(pars, 'gamma') && ~isempty(pars.gamma), pars.dw = pars.gamma; end
  if isfield(pars, 'u2')    && ~isempty(pars.u2),    pars.dw = pars.u2; end
  
  % constants
  SE2V = 437.393377;        % Convert sqrt(E)[meV] to v[m/s]
  V2K  = 1.58825361e-3;     % Convert v[m/s] to k[1/AA]
  K2V  = 1/V2K;
  VS2E = 5.22703725e-6;     % Convert (v[m/s])**2 to E[meV]
  
  % search in the vDOS data set in case missing properties are stored there
  % Ki, T, mass
  if isempty(pars.ki)
    pars.ki = Sqw_getT(gw, {'wavevector','momentum','IncidentWavevector','Ki'});
  end
  if isempty(pars.lambda)
    pars.lambda = Sqw_getT(gw, {'wavelength','lambda'});
  end
  if isempty(pars.ei)
    pars.ei = Sqw_getT(gw, {'IncidentEnergy','Ei'});
  end
  
  if isempty(pars.ki) && isfield(pars, 'lambda') && ~isempty(pars.lambda) && pars.lambda > 0
    pars.ki=2*pi/pars.lambda;
  end
  if isempty(pars.ki) && isfield(pars, 'ei') && ~isempty(pars.ei) && pars.ei > 0
    pars.ki=0.695*sqrt(pars.ei);
  end
  if isempty(pars.ki)
    pars.ki = 2.662;
  end
  
  if isempty(pars.t) || pars.t<=0
    pars.t = Sqw_getT(gw);
  end
  if isempty(pars.m) || pars.m<=0
    pars.m = Sqw_getT(gw, {'Masses','Molar_mass','Mass','Weight'});
  end
  if isempty(pars.phi)
    pars.phi = Sqw_getT(gw, {'DetectionAngles','Angle'});
  end
  if isempty(pars.phi)
    pars.phi = [ 10 120 ];
    disp([ mfilename ': using detector angular range phi=' mat2str(pars.phi) ' [deg] incident neutron wavelength.' ]);
  end
  % fail when missing information
  if isempty(pars.m) || pars.m<=0    
    disp([ mfilename ': WARNING: Unspecified molar mass (m). Using m=12 [g/mol]' ]); 
    pars.m = 12;
  end
  if isempty(pars.t) || pars.t<=0    
    disp([ mfilename ': WARNING: Unspecified temperature (T). Using T=10 [K]' ]); 
    pars.t = 10;
  end
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
  gw2   = gw1; setaxis(gw2, 1, -getaxis(gw2,1));
  gw    = cat(1, gw1, gw2);
  hw    = getaxis(gw,1);                % energy axis in [meV]
  
  % we search the smallest non zero energy value and create a 'delta' function
  [dw, index] =  min(abs(hw(:)));
  delta = 0*hw; delta(index) = 1;
  delta = iData(hw, delta);             % make it an iData for easier handling

  % compute the Debye-Waller function W(q) = <u2>/2. Schober p 328 -------------
  kT = pars.t/meVtoK;                        % T [K] -> [meV] = 11.6045
  
  % f(0) is an inverse energy. integrate -inf:inf
  % f(0)=\int(dw g(w)/w [n(w)+1]) Eq (10.58)
  nw = 1./(exp(hw./kT) -1);             % Bose population factor
  % compute the T1 term
  T1 = gw./hw.*(nw+1);                  % Eq. (10.68) p 329
  f0 = trapz(T1);
  if ~isempty(pars.dw) && pars.dw > 0
    % we use a specified DW 'gamma' value
    T1 = T1/f0;
    f0 = pars.dw/(q2toE/4/pars.m)/2;  % new DW factor sets f0=trapz(T1)
    T1 = T1*f0;
  end
  disp([ mfilename ': 1/f(0)=' num2str(1/f0) ' [meV]' ]);

  % Debye-Waller coefficient, aka gamma=2W/q^2 ; h/2m = q2toE/2/m
  W     = (q2toE/4/pars.m*f0);               % W = <u2> in Angs^2. integrate 0:inf
  gamma = 2*W;
  disp([ mfilename ': Debye-Waller coefficient gamma=<u^2>=' num2str(2*W) ...
    ' [Angs^2]. Debye-Waller function is 2W(q)=gamma*q^2. Debye-Waller factor is exp(-gamma.q^2)' ]);
  
  % compute the T[p] terms iteratively -----------------------------------------
  Tall = T1;
  
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
  
  % compute the symmetrized multi-phonon Tp terms ------------------------------
  Tsym = [];
  for p=1:pars.n
    if pars.n==1, Tp = Tall; else Tp = Tall(p); end
    Tpsym = exp(-hw./2/kT).*Tp;
    Tsym = [ Tsym Tpsym ];
  end
  
  % compute the integrated intensity -------------------------------------------
  Gw = [];
  % compute the dynamic range
  Qm   = [ Q(pars.ki, hw, min(pars.phi)) Q(pars.ki, hw, max(pars.phi)) ]; % size=[ numel(hw) 2 ]
  Qmax = max(Qm,[], 2);                               % size=[ numel(hw) 1 ]
  Qmin = min(Qm,[], 2);
  
  fact = 1;
  % compute the gDOS terms Eq (10.93)
  for p=1:pars.n
    if pars.n==1, Tpsym = Tsym; else Tpsym = Tsym(p); end
    fact = fact*p;  % !p
    dGw  = ((max(abs(pars.phi))-min(abs(pars.phi)))*1/8/pi/pars.ki.^2).*exp(hw./2/kT) ...
        .*(1/pars.m)^p/gamma^(p+1) ...
        .*(Ip(gamma*Qmax,p) - Ip(gamma*Qmin,p)).*Tpsym/fact;
    setalias(dGw, 'Temperature', pars.t, 'Temperature [K]');
    setalias(dGw, 'Weight',      pars.m, 'Molar masses [g/mol]');
    setalias(dGw, 'IncidentWavevector', pars.ki, 'neutron incident wavevector [Angs-1]');
    setalias(dGw, 'Wavelength', 2*pi/pars.ki, 'neutron incident wavelength [Angs]');
    setalias(dGw, 'DetectionAngles', pars.phi, 'Detection Angles [deg]');
    setalias(dGw, 'DebyeWallerCoeficient', gamma, ...
      'Debye-Waller coefficient gamma=<u^2> [Angs^2]. Debye-Waller function is 2W(q)=gamma*q^2. Debye-Waller factor is exp(-gamma.q^2)');
    dGw.Title= [ 'gDOS [p=' num2str(p) ']' ];
    title(dGw, [ 'gDOS [p=' num2str(p) '] from ' titl ]);
    dGw.Label= [ 'gDOS [p=' num2str(p) ']' ];
    dGw  = iData_vDOS(dGw);
    Gw   = [ Gw dGw ];
  end
  
  if nargout == 0
    fig=figure; 
    G1 = plus(Gw);        title(G1, 'gDOS [single+multi]'); G1.Label=title(G1); G1.Title=title(G1);
    G2 = plus(Gw(2:end)); title(G2, 'gDOS [multi]');        G2.Label=title(G2); G2.Title=title(G2);
    h=subplot([ G1 G2 ]); 
    set(fig, 'NextPlot','new');
  end

% ------------------------------------------------------------------------------
function ip = Ip(y, n)
  % Ip(y,p): compute Ip Q integrated terms
  % Ip(y) = int[ x^p exp(-x) dx ] from 0 to y
  %       = -y^p exp(-y) + p*I[p-1]
  %
  % input:
  %   y: value where to evaluate the function
  %   p: order to evaluate (recursively)
  
  % WARNING: there is an error in Schober Eq (10.83) and Reichard MUPHOCOR p3
  % Obviously, the proposed In formulae are all negative, but should be always 
  %   positive.
  % The recursion formula is wrong. The second term is +p*I[p-1]
  

  Ipm1= 1-exp(-y); % I0
  if n == 0, ip = Ipm1; return; end
  for p=1:n
    ip   = -y.^p.*exp(-y) + p*Ipm1;
    Ipm1 = ip;
  end
  ip(~isreal(ip)) = 0;
  
function q = Q(Ki, hw, phi)
  % Q: compute the Q value coresponding with a given scattering angle
  %
  % input:
  %   Ki:   incident neutron wavevector [Angs-1]
  %   w:    energy transfer [meV]
  %   phi:  scattering angle [deg]
  
  % constants
  SE2V = 437.393377;        % Convert sqrt(E)[meV] to v[m/s]
  V2K  = 1.58825361e-3;     % Convert v[m/s] to k[1/AA]
  K2V  = 1/V2K;
  VS2E = 5.22703725e-6;     % Convert (v[m/s])**2 to E[meV]
  Vi   = K2V*Ki;
  Ei   = VS2E*Vi.^2;
  % lambda = 2*pi./Ki;

  % compute Kf,Ef,Vf,Vi,Ei for the dynamic range computation
  Ef = Ei-hw;
  Ef(Ef<0) = 0; % neutron energy can not be negative
  Vf = SE2V*sqrt(Ef);
  Kf = V2K*Vf;
  phi= phi*pi/180;  % deg -> rad

  q2 = (Ki.^2 + Kf.^2 - cos(phi)*(2*Ki.*Kf));
  q  = sqrt(q2);
  q(Ef<=0) = 0;
  q  = q(:); % as a column with numel(hw)

 
  
% ------------------------------------------------------------------------------
  
% Impulse approx. Schober p 317 
% S(q,w) eq (10.25) Gaussian centred on recoil for Q -> Inf
% W(q)   eq (10.6)  for one single oscillator, isotropic scatterer
%        eq (10.3)  for a xtal

% sin(theta) integrated density of states. Oskotskii/Bredov eq (9.299) p 315
% incoherent approximation for coherent = 20% eq (9.295) p314
%   for an isotropic harmonic oscillator
