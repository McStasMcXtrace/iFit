function s = Sqw_phi2q(s)
% convert S(phi,w) to S(q,w). Requires wavelength

  if isempty(s), return; end
  [s,lambda] = Sqw_search_lambda(s);
  
  disp([ mfilename ': ' s.Tag ' ' s.Title ' Converting Axis 2 "' label(s, 2) '": angle [deg] to wavevector [Angs-1].' ]);
  Ei    = 81.805/lambda^2;
  phi   = s{2}; % angle
  hw    = s{1};
  if isvector(hw) && isvector(phi)
    s = meshgrid(s);
    phi   = s{2}; % angle
    hw    = s{1};
  end
  % we use: cos(phi) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf);
  Ei = 81.805/lambda^2; Ki = 2*pi./lambda; 
  Ef = Ei - hw;         Kf = sqrt(Ef/2.0721);
  q  = sqrt(Ki.^2 + Kf.^2 - 2*cos(phi*pi/180).*Ki.*Kf);

  s = setalias(s, 'q', q, 'Wavevector [Angs-1]');
  s = setaxis(s, 2, 'q');
