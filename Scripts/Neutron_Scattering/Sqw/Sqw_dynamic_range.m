function [s, sphi] = Sqw_dynamic_range(s, Ei, angles, options)
% [sqw_Ei,sphiw]=Sqw_dynamic_range(s,Ei): crop the S(|q|,w) to the available dynamic range
%   for given incident neutron energy.
%
%  The S(q,w) is a dynamic structure factor aka scattering function.
%
% The dynamic range is defined from the momwntum and energy conservation laws:
%  Ef         = Ei - w                                is positive
%  cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)  is within [-1:1]
%
% The incident neutron energy can be computed using:
%  Ei = 2.0721*Ki^2 = 81.8042/lambda^2 with Ki in [Angs-1] and lambda in [Angs]
%
% conventions:
% omega = Ei-Ef = energy lost by the neutron
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
%
% The scattering angle phi can be restricted to match a detection area
% with the syntax:
%   sqw_Ei=Sqw_dynamic_range(s, Ei, [angles])
%
% The syntax:
%   [sqw_Ei, sphiw] = Sqw_dynamic_range(s,...)
% also returns the S(phi,w) data set, which shows the detected signal vs scattering 
% angle and energy transfer.
%
% input:
%   s:  Sqw data set, e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
%   Ei: incoming neutron energy [meV]
%   angles: detection range in [deg] as a vector. Min and Max values are used.
% output:
%   sqw:   S(q,w)   cropped to dynamic range for incident energy Ei.
%   sphiw: S(phi,w) angular dynamic structure factor for incident energy Ei.
%
% Example: Sqw_dynamic_range(s, 14.8, [-20 135])
%
% See also: Sqw_Bosify, Sqw_deBosify, Sqw_symmetrize, Sqw_scatt_xs
% (c) E.Farhi, ILL. License: EUPL.

  if nargin == 0, return; end
  if ~isa(s, 'iData'), s=iData(s); end

  % first look at the given arguments
  if nargin < 2
    Ei = [];
  end
  if nargin < 3
    angles = [];
  end
  if nargin < 4, options=''; end
  sqw = [];

  % handle array of objects
  if numel(s) > 1
    for index=1:numel(s)
      sqw = [ sqw feval(mfilename, s(index), Ei, angles, options) ];
      s(index)=iData; % clear memory
    end
    s = sqw;
    return
  end
  
  % check if the energy range is limited
  w = s{1};
  if all(w(:) >= 0) || all(w(:) <= 0)
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ])
    disp([ '    seems to have its energy range w=[' num2str([ min(w(:)) max(w(:)) ]) '] defined only on one side.' ])
    disp('    Perhaps you should apply Sqw_symmetrize and Sqw_Bosify first ?');
  end
  
  % compute Ei, Ef, Ki, Kf
  % constants
  SE2V = 437.393377;        % Convert sqrt(E)[meV] to v[m/s]
  V2K  = 1.58825361e-3;     % Convert v[m/s] to k[1/AA]
  K2V  = 1/V2K;
  VS2E = 5.22703725e-6;     % Convert (v[m/s])**2 to E[meV]

  % compute Ki with given arguments
  if ~isempty(Ei)
    Ki = SE2V*V2K*sqrt(Ei);
  else
    % get dynamic limits from the given S(q,w)
    q = s{2};
    Ki = max(q(:))/2;
    Ei = max(abs(w(:)));
    if SE2V*V2K*sqrt(Ei) < Ki % check which of q,w limits the range
      Ki = SE2V*V2K*sqrt(Ei);
    end
  end
  
  if ~strcmp(options,'checked')
    s = Sqw_check(s); % in private
    if isempty(s), return; end
    
    % get the S(q,w) on a meshgrid (and possibly rebin/interpolate)
    s = meshgrid(s);
  end
  
  if numel(Ei) > 1
    for index=1:numel(Ei)
      sqw = [ sqw feval(mfilename, s, Ei(index), angles, 'checked') ];
    end
    s = sqw;
    return
  end
  
  % get axes
  w = s{1};
  q = s{2};

  % compute Vi
  if isempty(Ki), return; end
  Vi = K2V*Ki;
  % Ei = VS2E*Vi.^2;

  % compute Kf,Ef,Vf,Vi,Ei,lambda for the dynamic range computation
  lambda = 2*pi./Ki;

  Ef = Ei-w;
  
  Vf = SE2V*sqrt(Ef);
  Kf = V2K*Vf;

  fprintf(1, '%s: incoming: Ei=%g [meV] Ki=%g [Angs-1] lambda=%g [Angs] Vi=%g [m/s]\n', ...
    mfilename, Ei, Ki, lambda, Vi);
  fprintf(1, '  q=[%g:%g] w=[%g:%g]\n', min(q(:)), max(q(:)), min(w(:)), max(w(:)));

  % find: cos theta = (ki2 + kf2 - q2)/(2ki kf) is NOT between -1 and 1. theta=diffusion angle
  costheta = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf);
  if isempty(angles)
    index= find(abs(costheta) > 1 | Ef <= 0);
  else
    if isscalar(angles), angles = [ angles-1 angles+1 ]; end
    cost = cosd(angles);
    index= find(costheta < min(cost(:)) | costheta > max(cost(:)) | Ef <= 0);
  end
  s(index)=0;
  
  if isempty(angles)
    tt = sprintf(' Ei=%g meV', Ei);
  else
    tt = sprintf(' Ei=%g meV [%g:%g deg]', Ei, min(angles), max(angles));
  end
  
  s.DisplayName = strtrim([ s.DisplayName tt ]);
  s = setalias(s, 'IncidentEnergy', Ei, 'Incident Energy [meV]');
  if ~isempty(angles)
    s = setalias(s, 'DetectionAngles', angles, 'Detection Angles [deg]');
  end
  % compute the angles from 'q' values, and create an Angle alias
  if nargout > 1
    angles = real(acos(costheta))*180/pi; % this is e.g. a 2D array
    sphi = setalias(s, 'Angle', angles, 'Radial Angle [deg]');
    sphi = setaxis(sphi, 2, 'Angle');
  end
  
  s = title(s, strtrim([ s.Title tt ]));
  s.Label = strtrim(tt);
  
