function spt = qw2phit(self, varargin)
  % iFunc_Sqw2D: qw2phit: convert a S(q,w) into a S(phi,tof) Model (scattering angle)
  %
  % convert: iFunc_Sqw2D S(q,w) -> S(phi,t)
  %
  % An initial S(q,w) 2D Model is converted into a 2D (angle,tof) data set.
  %   The incident neutron energy is assumed to be monochromatic.
  %   The detector is assumed to be radially arranged (constant radius).
  %
  % New model parameters:
  %       Ei     Incident neutron energy [meV]
  %       Distance Sample-Detector radius [m]
  %
  %  The incident neutron energy can be computed using:
  %   Ei = 2.0721*Ki^2 = 81.8042/lambda^2 with Ki in [Angs-1] and lambda in [Angs]
  %
  % Example: s=sqw_visco_elastic_simple; spt=s.qw2phit; 
  %          plot(log(spt), [], 0:180, linspace(0, 0.003, 1000))
  %
  % spt = qw2phit(s)
  %
  % input:
  %   s:      Sqw 2D model with q as 1st axis (Angs-1), w as 2nd axis (meV).
  %
  % output:
  %   spt:    S(phi,w) 2D model [iFunc, phi in deg]
  
  index=numel(self.Parameters)+1;
  self.Parameters{end+1} = 'Ei Incident neutron energy [meV]';
  self.Parameters{end+1} = 'Distance Sample-Detector radius [m]';
  
  if isvector(self.Guess) && isnumeric(self.Guess)
    self.Guess = [ self.Guess(:)' 14.8 2.5 ];  % typical
  else
    if ~iscell(self.Guess), self.Guess = { self.Guess }; end
    self.Guess{end+1} = [ 14.8 2.5 ];
  end
  
  % the resulting 2D model will use phi[deg] and t[s]
  % prepend: from phi, we compute the q[Angs-1] axis
  prepend = { [ 'phi = x; t = y; Ei = p(' num2str(index) '); distance=p(' num2str(index+1) ');' ]; ...
    'Ki = sqrt(Ei/2.0721);' ;...
    'K2V=629.6224; Vi   = K2V*Ki;'; ...
    't_sample_detector = distance/Vi;'; ...
    'Ef = Ei.*(t_sample_detector./t).^2;'; ...
    'dtdE= t./(2.*Ef)*1e6'; ...
    'w  = Ei - Ef;';         ...
    'Kf = sqrt(Ef/2.0721);'; ...
    'q  = sqrt(Ki.^2 + Kf.^2 - 2*cosd(phi).*Ki.*Kf);'; ...
    'x  = q; y=w;'; ...
    
     };
  spt = prepend + self;
  % check for signal, and assign x axis as the initial angles.
  spt = spt + ...
    { 'signal(imag(signal) | imag(q))=0;'; ...
      'signal=signal.*dtdE; x=phi;y=t;' };
  spt.Name = [ mfilename '(' self.Name ')' ];
  spt = iFunc(spt); % check
  
end
