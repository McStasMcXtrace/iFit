function sqw = Sqw_phi2q(s, lambda, a_present, w_present)
% convert S(phi,w) to S(q,w). Requires wavelength

  sqw = [];
  if isempty(s), return; end
  if nargin < 2, lambda=[]; end
  if nargin < 3, a_present = 2; end
  if nargin < 4, w_present = 1; end
  
  if isempty(lambda)
    [s,lambda] = Sqw_search_lambda(s);
  end
  
  disp([ mfilename ': ' s.Tag ' ' s.Title ' Converting Axis ' num2str(a_present) ...
    ' "' label(s, a_present) '": angle [deg] to wavevector [Angs-1].' ]);
  Ei    = 81.805/lambda^2;
  phi   = getaxis(s,  a_present); % angle (assumed to be scattering angle)
  hw    = getaxis(s,  w_present);
  
  if isvector(hw) && isvector(phi)
    s = meshgrid(s);
    phi   = getaxis(s,a_present); % angle
    hw    = getaxis(s,w_present);
  end
  % we use: cos(phi) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf);
  Ei = 81.805/lambda^2; Ki = 2*pi./lambda; 
  Ef = Ei - hw;         
  if numel(Ef) < numel(phi)
    % Ef is a vector [Mx1], phi is a matrix [MxN]
    Ef = repmat(Ef, 1, size(phi,2));
  end
  Kf = sqrt(Ef/2.0721);
  q  = sqrt(Ki.^2 + Kf.^2 - 2*cos(phi*pi/180).*Ki.*Kf);

  sqw = copyobj(s);
  sqw = setalias(sqw, 'q', q, 'Wavevector [Angs-1]');
  sqw = setalias(sqw, 'IncidentWavelength', lambda);
  sqw = setaxis(sqw, a_present, 'q');
  
  sqw = commandhistory(sqw, 'phi2q', s, lambda);
  sqw.Label = 'S(q, w)';
  label(sqw, 0, [  'phi2q' '(' label(s, 0) ')' ]);
