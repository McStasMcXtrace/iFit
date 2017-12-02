function s = Sqw_q2phi(s)
% convert S(q,w) to S(phi,w). Requires wavelength

  if isempty(s), return; end
  
  s = iData(s);
  [s,lambda] = Sqw_search_lambda(s);
  
  disp([ mfilename ': ' s.Tag ' ' s.Title ' Converting Axis 2 "' label(s, 2) '": wavevector [Angs-1] to angle [deg].' ]);
  Ei    = 81.805/lambda^2;
  q     = s{2}; % wavevector
  hw    = s{1};
  if isvector(hw) && isvector(q)
    s = meshgrid(s);
    q     = s{2}; % angle
    hw    = s{1};
  end
  % we use: cos(phi) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf);
  Ei = 81.805/lambda^2; Ki = 2*pi./lambda; % direct geometry (monochromtic input)
  Ef = Ei - hw;         Kf = sqrt(Ef/2.0721);
  
  cos_phi = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf);
  phi     = acosd(cos_phi); % in degrees

  s = setalias(s, 'phi', phi, 'Scattering Angle [deg]');
  s = setaxis(s, 2, 'phi');
  
  
