function self=dynamic_range(self)
% iFunc_Sqw2D: dynamic_range: crop the S(|q|,w) to the available dynamic range
%    for given incident neutron energy.
% 
%  The dynamic range is defined from the momentum and energy conservation laws:
%   Ef         = Ei - w                                is positive
%   cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)  is within [-1:1]
%
% New model parameters:
%       Ei     Incident neutron energy [meV]
% 
%  The incident neutron energy can be computed using:
%   Ei = 2.0721*Ki^2 = 81.8042/lambda^2 with Ki in [Angs-1] and lambda in [Angs]
% 
%  conventions:
%  omega = Ei-Ef = energy lost by the neutron
%     omega > 0, neutron looses energy, can not be higher than p.ei (Stokes)
%     omega < 0, neutron gains energy, anti-Stokes
% 
%  syntax:
%    sqw_Ei        =dynamic_range(s)
% 
%  input:
%    s:      Sqw 2D model with q as 1st axis (Angs-1), w as 2nd axis (meV).
%  output:
%    sqw_Ei: S(q,w)   cropped to dynamic range for incident energy Ei.
% 
%  Example: s=sqw_visco_elastic; dynamic_range(s); plot(s);
% 
%  See also: iFunc_Sqw2D
%  (c) E.Farhi, ILL. license: EUPL.

  index=numel(self.Parameters)+1;
  self.Parameters{end+1} = 'Ei Incident neutron energy [meV]';
  self.Parameters{end+1} = 'Angle_min Minimum detector/scattering angle [deg]';
  self.Parameters{end+1} = 'Angle_max Maximum detector/scattering angle [deg]';
  
  if isvector(self.Guess) && isnumeric(self.Guess)
    self.Guess = [ self.Guess(:)' 14.8 5 135 ];  % typical
  else
    if ~iscell(self.Guess), self.Guess = { self.Guess }; end
    self.Guess{end+1} = [ 14.8 5 135 ];
  end
  t = self.Name;

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
  self.Name = [ mfilename '(' t ')' ];
  self = iFunc_Sqw2D(self); % check

  if nargout == 0 && ~isempty(inputname(1))
    assignin('caller',inputname(1),self);
  end

end % dynamic_range
