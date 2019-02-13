function spt = qw2phit(self, varargin)
  % iFunc_Sqw2D: qw2phit: convert a S(q,w) into a S(phi,tof) Model (scattering angle)
  %
  % convert: iFunc_Sqw2D S(q,w) -> S(phi,tof)
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
  %   spt:    S(phi,tof) 2D model [iFunc, phi in deg, tof in s]
  
  spt = copyobj(self);
  index=numel(spt.Parameters)+1;
  
  if isvector(spt.Guess) && isnumeric(spt.Guess)
    spt.Guess = [ spt.Guess(:)' 14.8 2.5 ];  % typical
  else
    if ~iscell(spt.Guess), spt.Guess = { spt.Guess }; end
    spt.Guess{end+1} = [ 14.8 2.5 ];
  end
  
  spt.Parameters{end+1} = [ 'Ei Incident neutron energy [meV] ' mfilename ];
  spt.Parameters{end+1} = [ 'Distance Sample-Detector radius [m] ' mfilename ];

  % the resulting 2D model will use phi[deg] and tof[s]
  % prepend: from phi, we compute the q[Angs-1] axis
  prepend = { [ 'phi_sav_phit = x; tof_sav_phit = y; Ei = p(' num2str(index) '); distance=p(' num2str(index+1) ');' ]; ...
    'Ki = sqrt(Ei/2.0721);' ;...
    'K2V=629.6224; Vi   = K2V*Ki;'; ...
    't_sample_detector = distance/Vi;'; ...
    'Ef = Ei.*(t_sample_detector./tof_sav_phit).^2;'; ...
    'dtdE= tof_sav_phit./(2.*Ef)*1e6'; ...
    'w  = Ei - Ef;';         ...
    'Kf = sqrt(Ef/2.0721);'; ...
    'q  = sqrt(Ki.^2 + Kf.^2 - 2*cosd(phi_sav_phit).*Ki.*Kf);'; ...
    'x  = q; y=w;'; ...
    
     };
  spt = prepend + spt;
  % check for signal, and assign x axis as the initial angles.
  spt = spt + ...
    { 'signal(imag(signal) | imag(q))=0;'; ...
      'signal=signal.*dtdE; x=phi_sav_phit; y=tof_sav_phit;' };
  spt.Name = [ mfilename '(' self.Name ')' ];
  spt = iFunc(spt); % check
  
end
