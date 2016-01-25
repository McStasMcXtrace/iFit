function sigma=Sqw_moments(data, classical)
% Sqw_moments: compute harmonic frequencies
%
% input:
%   data: iData object for S(q,w)
%   classical: 0 for non symmetric S(q,w) [with Bose], 1 for symmetric (from MD)

  data = Sqw_check(data);
  if isempty(data), sigma=[]; return; end

  kT      = Sqw_getT(data)/11.604;   % kbT in meV;
  M       = data.Data.weight;               % mass
  
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
  M2      = abs(trapz(w.^2.*data)); % M2 cl = wc^2
  M3      = abs(trapz(abs(w).^3.*data));
  M4      = abs(trapz(w.^4.*data));
  
  % half width from normalized 2nd frequency moment J-P.Hansen and I.R.McDonald 
  % Theory of simple liquids Academic Press New York 2006.
  wq      = 2*q.*sqrt(kT./M0/M);  % Lovesey p180 eq. 5.38 = w0
  
  if classical
    % all odd moments are 0, even are to be multiplied by 2 when using S(q,w>0)
    % M2 = q.^2.*kT/M
    wc      = sqrt(M2./M0); % sqrt(<w²S>/s(q)) == sqrt(q kT/M/s(q))
    wl      = M3./M2; % maxima wL(q) of the longitudinal current correlation function ~ wl
  else
    wc      = sqrt(2*kT.*M1./M0); 
    wl      = sqrt(M3./M1); 
  end
  wc.Label='wc collective dispersion';  wl.Label='wl harmonic excitation'; wq.Label='w_q=q \sqrt(k_B T/m S(q))';

  sigma.wc=wc;
  sigma.wl=wl;
  sigma.wq=wq;
  return
  % now fit gaussians for each q value...
  for index=1:length(q)
    this  = data(index,:);
    w.dwG(index) = std(this);
    p = [ max(this) 0 dwG(index) 0 ];
    p = fits(this, 'gauss',p,'',[ 0 1 0 0]);
    w.dwF(index) = p(3);
  end
