function self=dynamic_range(self)
% iFunc_Sqw2D: dynamic_range: crop the S(|q|,w) to the available dynamic range
%    for given incident neutron energy.
% 
%  The dynamic range is defined from the momentum and energy conservation laws:
%   Ef         = Ei - w                                is positive
%   cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)  is within [-1:1]
% 
%  The incident neutron energy can be computed using:
%   Ei = 2.0721*Ki^2 = 81.8042/lambda^2 with Ki in [Angs-1] and lambda in [Angs]
% 
%  conventions:
%  omega = Ei-Ef = energy lost by the neutron
%     omega > 0, neutron looses energy, can not be higher than p.ei (Stokes)
%     omega < 0, neutron gains energy, anti-Stokes
% 
%  The scattering angle phi can be restricted to match a detection area
%  with the syntax:
%    sqw_Ei=dynamic_range(s, Ei, [angles])
% 
%  The syntax:
%    [sqw_Ei, sphiw] = dynamic_range(s,...)
%  also returns the S(phi,w) data set, which shows the detected signal vs scattering 
%    angle and energy transfer.
% 
%   Input arguments can be given in order, or with name-value pairs, or as a 
%     structure with named fields. In addition to Ei, one can specify 'lambda' or 'Ki'
%   e.g.:
%    sqw_Ei=dynamic_range(s, 'lambda', 2.36);
%    sqw_Ei=dynamic_range(s, 'Ki', 2.665);
% 
%  syntax:
%    sqw_Ei        =dynamic_range(s, Ei)
%    [sqw_Ei,sphiw]=dynamic_range(s, Ei)
%    [sqw_Ei,sphiw]=dynamic_range(s, Ei, [min_angle max_angle])
%    [sqw_Ei,sphiw]=dynamic_range(s, 'Ei', Ei, 'angles', [min_angle max_angle])
%    [sqw_Ei,sphiw]=dynamic_range(s, 'lambda', lambda)
% 
%  input:
%    s:      Sqw 2D data set with q as 1st axis (Angs-1), w as 2nd axis (meV).
%    Ei:     incoming neutron energy [meV]
%    angles: scattering detection range in [deg] as a vector. Min and Max values are used.
%    'lambda','Ki': additional named arguments to specify the incident energy
%  output:
%    sqw_Ei: S(q,w)   cropped to dynamic range for incident energy Ei.
%    sphiw:  S(phi,w) angular dynamic structure factor for incident energy Ei.
% 
%  Example: s=sqw_visco_elastic; dynamic_range(s)
% 
%  See also: iFunc_Sqw2D
%  (c) E.Farhi, ILL. license: EUPL.

  index=numel(self.Parameters)+1;
  self.Parameters{end+1} = 'Ei Incident neutron energy [meV]';
  self.Parameters{end+1} = 'Angle_min Minimum detector/scattering angle [deg]';
  self.Parameters{end+1} = 'Angle_max Maximum detector/scattering angle [deg]';
  
  if isvector(self.Guess) && isnumeric(self.Guess)
    self.Guess = [ self.Guess(:)' 14.8 5 135 ];  % typical
  end

  self.Expression{end+1} = [ 'q = x; w = y; Ei = p(' num2str(index) ');' ];
  self.Expression{end+1} = [ 'angles = p([' num2str(index+1:index+2) ']);' ];  
  self.Expression{end+1} = 'if prod(angles) < 0, angles = [ 0 max(abs(angles)) ]; end';
  self.Expression{end+1} = 'SE2V = 437.393377;        % Convert sqrt(E)[meV] to v[m/s]';
  self.Expression{end+1} = 'V2K  = 1.58825361e-3;     % Convert v[m/s] to k[1/AA]';
  self.Expression{end+1} = 'Ki = SE2V*V2K*sqrt(Ei);';
  self.Expression{end+1} = 'Ef = Ei-w;';
  self.Expression{end+1} = 'Kf = SE2V*V2K*sqrt(Ef);';
  self.Expression{end+1} = 'costheta = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf);';
  self.Expression{end+1} = 'cost = cosd(angles);';
  self.Expression{end+1} = 'index= find(costheta < min(cost(:)) | costheta > max(cost(:)) | Ef <= 0);';
  self.Expression{end+1} = 'signal(index) = 0; % actually set the dynamic range';
  
  self = iFunc_Sqw2D(self); % check

  if nargout == 0 && ~isempty(inputname(1))
    assignin('caller',inputname(1),self);
  end

end % dynamic_range
