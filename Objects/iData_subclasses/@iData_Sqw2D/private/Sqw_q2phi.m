function spw = Sqw_q2phi(s, lambda)
% convert S(q,w) to S(phi,w). Requires wavelength

  spw = [];
  if isempty(s), return; end
  if nargin < 2, lambda=[]; end
  
  s = iData(s);
  if isempty(lambda)
    [s,lambda] = Sqw_search_lambda(s);
  end

  if isempty(lambda)
    disp([ mfilename ': ' s.Tag ' ' s.Title ' Using lambda=2.36 [Angs].' ]);
    lambda = 2.36;
  end
  
  disp([ mfilename ': ' s.Tag ' ' s.Title ' Converting Axis 2 "' label(s, 2) ...
    '": wavevector [Angs-1] to scattering angle [deg].' ]);
  Ei    = 81.805./lambda.^2;
  q     = s{2}; % wavevector
  hw    = s{1};
  if isvector(hw) && isvector(q)
    s = meshgrid(s);
    q     = s{2}; % angle
    hw    = s{1};
  end
  % we use: cos(phi) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf);
  Ki = 2*pi./lambda; % direct geometry (monochromatic input)
  Ef = Ei - hw;         Kf = sqrt(Ef/2.0721);

  cos_phi = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf);
  cos_phi(Ef < 0 | abs(cos_phi) > 1) = nan;
  phi     = acosd(cos_phi); % in degrees

  spw = iData(copyobj(s)); % make it a true iData
  spw = setalias(spw, 'phi', phi, 'Scattering Angle [deg]');
  spw = setalias(spw, 'IncidentWavelength', lambda);
  spw = setaxis(spw, 2, 'phi');
  
  spw = commandhistory(spw, 'qw2phiw', s, lambda);
  spw.Label = 'S(phi, w)';
  label(spw, 0, [  'qw2phiw [q2phi]' '(' label(s, 0) ')' ]);
  
  
