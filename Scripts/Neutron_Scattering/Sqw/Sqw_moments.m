function sigma=Sqw_moments(data, M, T, classical)
% moments=Sqw_moments(sqw, M, T, classical): compute Sqw moments (harmonic frequencies)
%
%   Compute the structure factor (moment 0), recoil energy (moment 1) and the
%     collective, harmonic and mean energy transfer dispersions.
%
% The result is given as an iData array with data sets:
%   S(q) = \int S(q,w) dw = <S(q,w)>                 structure factor [moment 0]
%   Er   = \int w*S(q,w) dw = <wS(q,w)> = h2q2/2M       recoil energy [moment 1]
%   Wc   = sqrt(2kT*Er/S(q))                    collective/isothermal dispersion
%   Wl                                          harmonic/longitudinal excitation
%   Wq   = 2q*sqrt(kT/S(q)/M)                               mean energy transfer
%   M2   = <w2S(q,w)>                                                 [moment 2]
%   M3   = <w3S(q,w)>                                                 [moment 3]
%   M4   = <w4S(q,w)>                                                 [moment 4]
%
% input:
%   data: Sqw data set e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
%   M: molar weight of the material atom/molecule in [g/mol].
%     when omitted or empty, it is searched 'weight' or 'mass' is the object.
%   T: when given, Temperature to use. When not given or empty, the Temperature
%      is searched in the object. The temperature is in [K]. 1 meV=11.605 K.
%   classical: 0 for non symmetric S(q,w) [with Bose, from exp.], 1 for symmetric (from MD)
%     when omitted or empty, this is guessed from the data set when possible
%
% output:
%   moments=[ sq M1 wc wl wq M2 M3 M4 ] as iData array
%
% Reference: 
%   Helmut Schober, Journal of Neutron Research 17 (2014) pp. 109
%   Lovesey, Theory of Neutron Scattering from Condensed Matter, Vol 1, p180 eq. 5.38 (w0)
%   J-P.Hansen and I.R.McDonald, Theory of simple liquids Academic Press New York 2006.
%
% Example: moments = sqw_moments(iData(fullfile(ifitpath,'Data','SQW_coh_lGe.nc')))
% (c) E.Farhi, ILL. License: EUPL.

  sigma = [];
  if nargin == 0, return; end
  if ~isa(data, 'iData'), data=iData(data); end
  if nargin < 2, M = []; end
  if nargin < 3, classical = []; end
  if nargin < 4, T = []; end
  
  data = Sqw_check(data);
  if isempty(data), sigma=[]; return; end
  
  % guess when omitted arguments
  if isempty(classical) && (isfield(data,'classical') || ~isempty(findfield(data, 'classical')))
    classical = data.classical;
  end
  if isempty(M) && isfield(data, 'weight')
    M       = data.weight;               % mass
  end
  if isempty(T)
    T = Sqw_getT(data);
  end
  
  % check input parameters
  if isempty(classical)
    disp([ mfilename ': ERROR: The data set ' data.Tag ' ' data.Title ' from ' data.Source ])
    disp('   does not provide information about classical/quantum data set.');
    disp('   Use Sqw_moments(data, M, T, classical=0 or 1)');
    return
  end
  
  if isempty(T) || T<=0
    disp([ mfilename ': ERROR: Temperature undefined: The data set ' data.Tag ' ' data.Title ' from ' data.Source ]);
    disp('    does not have any temperature defined. Use Sqw_moments(data, M, T, classical).' );
    return
  end
  
  kT      = T/11.604;   % kbT in meV;
  q       = data{2};
  w       = data{1}; 
  
  % clean low level data
  i=find(log(data) < -15);
  data(i) = 0;
  
  sq      = abs(trapz(data)); % S(q) from the data itself
  M0      = sq;
  % w2R = 2 kT M1
  % w2R 1/2/kT = wS = M1 and w0^2 = 1/S(q) w2R = 1/S(q) 2 kT M1 = q2 kT/M/M0
  M1      = abs(trapz(abs(w).*data));    % = h2q2/2/M recoil when non-classical, 0 for classical symmetrized
  if ~classical && isempty(M)
    % try to extract a mass from the recoil
    mn      = 1.675E-027;      % neutron mass [kg]
    e       = 1.602E-019;      % [C]
    HBAR    = 1.05457168e-34;  % Plank/2PI
    kb      = 1.381E-023;      % Boltzmann [J/K]
    q2toE   = HBAR*HBAR/2/mn/e*1000*1e20; % = 2.0723 = [Angs^-2] to [meV] 
    C       = e/1000/kb/T;
    M=mean(q.*q*2.0723./M1);
    if M >= 1 && M < 2000
      disp([ mfilename ': INFO: The data set ' data.Tag ' ' data.Title ' from ' data.Source  ]);
      disp([ '    Recoil provides a mass M=' num2str(M) ' [g/mol]. Wq may be wrong.' ]);
    else
      M = [];
    end
  end
  if isempty(M)
    disp([ mfilename ': WARNING: Mass undefined: The data set ' data.Tag ' ' data.Title ' from ' data.Source ]);
    disp('    does not provide information about the material molar weight. Use Sqw_moments(data, M).')
    disp('    Ignoring: The Wq frequency will be empty.');
  end
  M2      = abs(trapz(w.^2.*data)); % M2 cl = wc^2
  M3      = abs(trapz(abs(w).^3.*data));
  M4      = abs(trapz(w.^4.*data));
  
  % half width from normalized 2nd frequency moment J-P.Hansen and I.R.McDonald 
  % Theory of simple liquids Academic Press New York 2006.
  if ~isempty(M) && isnumeric(M)
    wq      = 2*q.*sqrt(kT./M0/M);  % Lovesey p180 eq. 5.38 = w0
    wq.Label='w_q=q \surd kB T/m S(q) mean energy transfer';
  else
    wq = iData;
  end
  
  if classical
    % all odd moments are 0, even are to be multiplied by 2 when using S(q,w>0)
    % M2 = q.^2.*kT/M
    wc      = sqrt(M2./M0); % sqrt(<w2S>/s(q)) == q sqrt(kT/M/s(q)) collective/isothermal
    wl      = M3./M2; % maxima wL(q) of the longitudinal current correlation function ~ wl
  else
    wc      = sqrt(2*kT.*M1./M0); 
    wl      = sqrt(M3./M1); 
  end
  
  sq.Label='S(q) structure factor';
  M1.Label='recoil E_r=h^2q^2/2M <wS> 1st moment';
  wc.Label='w_c collective/isothermal dispersion';  
  wl.Label='w_l harmonic/longitudinal excitation';
  M2.Label='<w2S> 2nd moment';
  M3.Label='<w3S> 3rd moment';
  M4.Label='<w4S> 4th moment';
  
  sigma =[ sq M1 wc wl wq M2 M3 M4 ];

  return
  
  
  % now fit gaussians for each q value...
  for index=1:length(q)
    this  = data(index,:);
    w.dwG(index) = std(this);
    p = [ max(this) 0 dwG(index) 0 ];
    p = fits(this, 'gauss',p,'',[ 0 1 0 0]);
    w.dwF(index) = p(3);
  end
