function spw = qw2phiw(self, varargin)
  % iFunc_Sqw2D: qw2phiw: convert a S(q,w) into a S(phi,w) Model (scattering angle)
  %
  % convert: iFunc_Sqw2D S(q,w) -> S(phi,w)
  %
  % An initial S(q,w) 2D Model is converted into a 2D (angle,w) data set.
  %   The incident neutron energy is assumed to be monochromatic.
  %
  % New model parameters:
  %       Ei     Incident neutron energy [meV]
  %
  %  The incident neutron energy can be computed using:
  %   Ei = 2.0721*Ki^2 = 81.8042/lambda^2 with Ki in [Angs-1] and lambda in [Angs]
  %
  % Example: s=sqw_visco_elastic_simple; spw=s.qw2phiw; 
  %          plot(log(spw), [], 0:180, -30:30)
  %
  % spw = qw2phiw(s)
  %
  % input:
  %   s:      Sqw 2D model with q as 1st axis (Angs-1), w as 2nd axis (meV).
  %
  % output:
  %   spw:    S(phi,w) 2D model [iFunc, phi in deg]
  
  spw = copyobj(self);
  index=numel(spw.Parameters)+1;
  spw.Parameters{end+1} = [ 'Ei Incident neutron energy [meV] ' mfilename ];
  
  if isvector(spw.Guess) && isnumeric(spw.Guess)
    spw.Guess = [ spw.Guess(:)' 14.8 ];  % typical
  else
    if ~iscell(spw.Guess), spw.Guess = { spw.Guess }; end
    spw.Guess{end+1} = [ 14.8 ];
  end
  
  % the resulting 2D model will use phi[deg] and w[meV]
  % prepend: from phi, we compute the q[Angs-1] axis
  prepend = { [ 'phi_sav_phit = x; w = y; Ei = p(' num2str(index) ');' ]; ...
    'Ef = Ei - w;';          ...
    'Ki = sqrt(Ei/2.0721);' ;...
    'Kf = sqrt(Ef/2.0721);'; ...
    'q  = sqrt(Ki.^2 + Kf.^2 - 2*cosd(phi_sav_phit).*Ki.*Kf);'; ...
    'x  = q;' };
  spw = prepend + self;
  % check for signal, and assign x axis as the initial angles.
  spw = spw + ...
    { 'signal(imag(q) | imag(signal))=0;'; ...
      'x=phi_sav_phit;' };
  spw.Name = [ mfilename '(' self.Name ')' ];
  spw = iFunc(spw); % check
  
end
