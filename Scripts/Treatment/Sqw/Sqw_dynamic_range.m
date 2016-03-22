function s = Sqw_dynamic_range(s, Ei, angles, options)
% Sqw_dynamic_range(s,Ei): crop the S(|q|,w) to the available dynamic range
%   for given incident neutron energy.
%
%  The S(q,w) is a dynamic structure factor aka scattering function.
%
% The dynamic range is defined from the momwntum and energy conservation laws:
%  Ef         = Ei - w                                is positive
%  cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)  is within [-1:1]
%
% The positive energy values in the S(q,w) map correspond to Stokes processes, 
% i.e. material gains energy, and neutrons loose energy when scattered.
%
% The scattering angle theta can be restricted to match a detection area
% with the syntax:
%   Sqw_dynamic_range(s, Ei, [angles])
%
% input:
%   s:  Sqw data set, with q as 1st axis, w as 2nd axis
%   Ei: incoming neutron energy [meV]
%   angles: detection range in [deg] as a vector. Min and Max values are used.
% output:
%   sqw: S(q,w) cropped to dynamic range
%
% Example: Sqw_dynamic_range(s, 14.8, [-20 135])
%
% See also: Sqw_Bosify, Sqw_deBosify, Sqw_symmetrize, Sqw_scatt_xs

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
  if any(Ef < 0)
      % the definition of the energy is reverted: Ef = Ei+w
      Ef = Ei+w;
  end
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
    cost = cosd(angles);
    index= find(costheta < min(cost(:)) | costheta > max(cost(:)) | Ef <= 0);
  end
  s(index)=0;
  s.DisplayName = strtrim([ s.DisplayName ' Ei=' num2str(Ei) ]);
  s = title(s, strtrim([ s.Title ' Ei=' num2str(Ei) ]));
  
