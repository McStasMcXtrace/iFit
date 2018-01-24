function [Gw, Tsym] = multi_phonons(gw, Ki, T, sigma, m, phi, n)
% iData_vDOS: multi_phonons: compute the integrated multi-phonon DOS from an initial density of states
%
% The 'generalized' neutron weighted density of states (gDOS) is computed from an 
% initial vibrational density of states (vDOS).
% Missing arguments (or given as [] empty), are searched in the initial density 
% of states.
%
% This implementation is in principle exact for an isotropic monoatomic material,
% e.g. a liquid or powder.
% This methodology is a complete rewrite of the MUPHOCOR code.
%
% Reference:
%   H. Schober, Journal of Neutron Research 17 (2014) 109–357
%     DOI 10.3233/JNR-140016 (see esp. pages 328-331)
%   V.S. Oskotskii, Sov. Phys. Solid State 9 (1967), 420.
%   A. Sjolander, Arkiv for Fysik 14 (1958), 315.
%   W. Reichardt, MUPHOCOR Karlsruhe Report 13.03.01p06L (1984)
%
% syntax:
%   Gw = multi_phonons(gw, Ki, T, sigma, m, n)
%
% input:
%   gw:   the vibrational density of states per [meV] [iData]
%   Ki:   incident wavevector [Angs-1]. Ki=2*pi/lambda=0.695*sqrt(Ei)
%   T:    temperature [K]
%   sigma: neutron cross section [barns]
%   m:    mass [g/mol]
%   phi:  detector angles, min and max [deg]
%   n:    number of iterations in the series expansion, e.g. 5
%
% output:
%   Gw:   neutron weighted gDOS terms, to be summed [iData_vDOS array]
%   Wq:   half Debye-Waller factor. The DW function is exp(-2*Wq)  [iData vs q]
%   Tp:   p-phonon terms [iData array]
%
% Example: s=sqw_cubic_monoatomic; multi_phonons(dos(s))
%
% See also: iData_Sqw2D/dos, iData_vDOS/multi_phonons_incoherent
% (c) E.Farhi, ILL. License: EUPL.

  % check input parameters
  if nargin < 2, Ki   =[]; end
  if nargin < 3, T    =[]; end
  if nargin < 4, sigma=[]; end
  if nargin < 5, m    =[]; end
  if nargin < 6, phi  =[]; end
  if nargin < 7, n    =[]; end
  
  % constants
  SE2V = 437.393377;        % Convert sqrt(E)[meV] to v[m/s]
  V2K  = 1.58825361e-3;     % Convert v[m/s] to k[1/AA]
  K2V  = 1/V2K;
  VS2E = 5.22703725e-6;     % Convert (v[m/s])**2 to E[meV]
  
  % search in the vDOS data set in case missing properties are stored there
  % Ki, T, mass, sigma
  if isempty(Ki)
    Ki = Sqw_getT(gw, {'wavevector','momentum','IncidentWavevector','Ki'});
    if isempty(Ki) || Ki<=0
      lambda = Sqw_getT(gw, {'wavelength','lambda'});
      if ~isempty(lambda) && lambda>0
        Ki = 2*pi/lambda;
      else
        Ei = Sqw_getT(gw, {'IncidentEnergy','Ei'});
        if isempty(Ei) || Ei<=0
          w  = getaxis(gw, 1);
          Ei = max(abs(w(:)));
          disp([ mfilename ': using Ei=' num2str(Ei) ' [meV] incident neutron energy.' ]);
        end
        Ki = SE2V*V2K*sqrt(Ei);
      end
    end
  end
  if isempty(T) || T<=0
    T = Sqw_getT(gw);
  end
  if isempty(m) || m<=0
    m = Sqw_getT(gw, {'Masses','Molar_mass','Mass','Weight'});
  end
  if isempty(sigma) || all(sigma<=0)
    sigma = Sqw_getT(gw, {'sigma_coh','sigma_inc','sigma'});
  end
  if isempty(phi)
    phi = Sqw_getT(gw, {'DetectionAngles','Angle'});
  end
  if isempty(phi)
    phi = [ 10 120 ];
    disp([ mfilename ': using detector angular range phi=' mat2str(phi) ' [deg] incident neutron wavelength.' ]);
  end
  % fail when missing information
  if isempty(m) || m<=0    error([ mfilename ': Unspecified molar mass (m). Use multi_phonons(g, Ki, T, sigma, m)' ]); end
  if isempty(T) || T<=0    error([ mfilename ': Unspecified temperature (T). Use multi_phonons(g, Ki, T)' ]); end
  if isempty(sigma) disp([ mfilename ': WARNING: Unspecified scattering cross section (sigma). Using sigma=1 barn. Scale result accordingly.' ]); sigma = 1; end
  if isempty(n),    n=5; end
  
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
  kT = T/meVtoK;                        % T [K] -> [meV] = 11.6045
  
  % f(0) is an inverse energy. integrate -inf:inf
  % f(0)=\int(dw g(w)/w [n(w)+1]) Eq (10.58)
  nw = 1./(exp(hw./kT) -1);             % Bose population factor
  % compute the T1 term
  T1 = gw./hw.*(nw+1);                  % Eq. (10.68) p 329
  f0 = trapz(T1);
  disp([ mfilename ': 1/f(0)=' num2str(1/f0) ' [meV]' ]);

  % Debye-Waller coefficient, aka gamma=2W/q^2 ; h/2m = q2toE/2/m
  W     = (q2toE/4/m*f0);               % W = <u2> in Angs^2. integrate 0:inf
  gamma = 2*W;
  disp([ mfilename ': Debye-Waller coefficient gamma      =' num2str(2*W) ...
    ' [Angs^2]. Debye-Waller function is 2W(q)=gamma*q^2' ]);
  
  % compute the T[p] terms iteratively -----------------------------------------
  Tall = T1;
  
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
  
  % compute the symmetrized multi-phonon Tp terms ------------------------------
  Tsym = [];
  for p=1:n
    if n==1, Tp = Tall; else Tp = Tall(p); end
    Tpsym = exp(-hw./2/kT).*Tp;
    Tsym = [ Tsym Tpsym ];
  end
  
  % compute the integrated intensity -------------------------------------------
  Gw = [];
  % compute the dynamic range
  Qm   = [ Q(Ki, hw, min(phi)) Q(Ki, hw, max(phi)) ]; % size=[ numel(hw) 2 ]
  Qmax = max(Qm,[], 2);                               % size=[ numel(hw) 1 ]
  Qmin = min(Qm,[], 2);
  
  fact = 1;
  % compute the gDOS terms Eq (10.93)
  for p=1:n
    if n==1, Tpsym = Tsym; else Tpsym = Tsym(p); end
    fact = fact*p;  % !p
    dGw  = ((max(abs(phi))-min(abs(phi)))*sigma/8/pi/Ki.^2).*exp(hw./2/kT) ...
        .*(1/m)^p/gamma^(p+1) ...
        .*(Ip(gamma*Qmax,p) - Ip(gamma*Qmin,p)).*Tpsym/fact;
    setalias(dGw, 'Temperature', T, 'Temperature [K]');
    setalias(dGw, 'Weight',      m, 'Molar masses [g/mol]');
    setalias(dGw, 'Sigma',   sigma, 'Neutron scattering cross section [barns]');
    setalias(dGw, 'IncidentWavevector', Ki, 'neutron incident wavevector [Angs-1]');
    setalias(dGw, 'Wavelength', 2*pi/Ki, 'neutron incident wavelength [Angs]');
    setalias(dGw, 'DetectionAngles', phi, 'Detection Angles [deg]');
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
