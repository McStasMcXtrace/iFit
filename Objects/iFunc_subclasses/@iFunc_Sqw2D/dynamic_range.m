function self=dynamic_range(self, method)
% iFunc_Sqw2D: dynamic_range: crop the S(|q|,w) to the available dynamic range
%    for given incident neutron energy.
% 
%  The dynamic range is defined from the momentum and energy conservation laws:
%   Ef         = Ei - w                                is positive
%   cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)  is within [-1:1]
%
% The dynamic range can be computed in direct (constant Ei) or indirect
% (constant Ef) mode. Default is 'direct'.
%
% New model parameters:
%       Ei     Incident neutron energy [meV] in direct mode
% or
%       Ef     Incident neutron energy [meV] in indirect mode
% 
%  The incident neutron energy can be computed using e.g.:
%   Ei = 2.0721*Ki^2 = 81.8042/lambda^2 with Ki in [Angs-1] and lambda in [Angs]
% 
%  conventions:
%  omega = Ei-Ef = energy lost by the neutron
%     omega > 0, neutron looses energy, can not be higher than p.ei (Stokes)
%     omega < 0, neutron gains energy, anti-Stokes
% 
%  syntax:
%    sqw_Ei        =dynamic_range(s)
%    sqw_Ef        =dynamic_range(s, 'indirect')
% 
%  input:
%    s:      Sqw 2D model with q as 1st axis (Angs-1), w as 2nd axis (meV).
%  output:
%    sqw_Ei: S(q,w)   cropped to dynamic range for incident energy Ei.
% 
%  Example: s=sqw_visco_elastic; dynamic_range(s); plot(s);
% 
%  See also: iFunc_Sqw2D
%  

  if nargin < 2, method = ''; end
  if isempty(method), method='direct'; end

  % add new Parameters
  if isvector(self.Guess) && isnumeric(self.Guess)
    self.Guess = [ self.Guess(:)' 14.8 5 135 ];  % typical
  else
    if ~iscell(self.Guess), self.Guess = { self.Guess }; end
    self.Guess{end+1} = [ 14.8 5 135 ];
  end
  
  switch lower(method)
  case {'indirect','ef'}
  self.Parameters{end+1} = [ 'Ef Final neutron energy [meV] ' mfilename ];
  otherwise % {'direct','default','defaults', 'ei'}
  self.Parameters{end+1} = [ 'Ei Incident neutron energy [meV] ' mfilename ];
  end
  self.Parameters{end+1} = [ 'Angle_min Minimum detector/scattering angle [deg] ' mfilename ];
  self.Parameters{end+1} = [ 'Angle_max Maximum detector/scattering angle [deg] ' mfilename ];

  % assemble expression
  t = self.Name;
  index=numel(self.Parameters);
  self.Expression{end+1} = [ 'q = x; w = y;' ];
  self.Expression{end+1} = [ 'angles = p([' num2str(index-1:index) ']);' ];  
  self.Expression{end+1} = 'if prod(angles) < 0, angles = [ 0 max(abs(angles)) ]; end';
  self.Expression{end+1} = 'SE2V = 437.393377;        % Convert sqrt(E)[meV] to v[m/s]';
  self.Expression{end+1} = 'V2K  = 1.58825361e-3;     % Convert v[m/s] to k[1/AA]';
  switch lower(method)
  case {'indirect','ef'}
  self.Expression{end+1} = ['Ef = p(' num2str(index-2) ');' ];
  self.Expression{end+1} = 'Kf = SE2V*V2K*sqrt(Ef);';
  self.Expression{end+1} = 'Ei = Ef+w;';
  self.Expression{end+1} = 'Ki = SE2V*V2K*sqrt(Ei);';
  self.Expression{end+1} = 'this.UserData.Kf=Kf; this.UserData.lambda=2*pi/Kf;';
  otherwise % {'direct','default','defaults', 'ei'}
  self.Expression{end+1} = ['Ei = p(' num2str(index-2) ');' ];
  self.Expression{end+1} = 'Ki = SE2V*V2K*sqrt(Ei);';
  self.Expression{end+1} = 'Ef = Ei-w;';
  self.Expression{end+1} = 'Kf = SE2V*V2K*sqrt(Ef);';
  self.Expression{end+1} = 'this.UserData.Ki=Ki; this.UserData.lambda=2*pi/Ki;'
  end
  self.Expression{end+1} = 'costheta = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf);';
  self.Expression{end+1} = 'cost = cosd(angles);';
  self.Expression{end+1} = 'index= find(costheta < min(cost(:)) | costheta > max(cost(:)) | Ef <= 0 | Ei <= 0);';
  self.Expression{end+1} = 'signal(index) = 0; % actually set the dynamic range';
  self.Name = [ mfilename '(' t ')' ];
  self.UserData.dynamic_range_mode = method;
  self = iFunc_Sqw2D(self); % check

  if nargout == 0 && ~isempty(inputname(1))
    assignin('caller',inputname(1),self);
  end

end % dynamic_range
