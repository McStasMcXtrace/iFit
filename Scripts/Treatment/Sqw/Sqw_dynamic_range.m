function sqw = Sqw_dynamic_range(s, Ei, angles)
% Sqw_dynamic_range(s,Ei): crop the S(|q|,w) to the available dynamic range
%   for given incident neutron energy.
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
% See also: Sqw_Bosify, Sqw_deBosify, Sqw_symmetrize, Sqw_total_xs

  % first look at the given arguments
  if nargin < 2
    Ei = [];
  end
  if nargin < 3
    angles = [];
  end

  % handle array of objects
  if numel(s) > 1
    sqw = [];
    for index=1:numel(s)
      sqw = [ sqw feval(mfilename, s(index), Ei, angles) ];
    end
    return
  end

  % check if the data set is Sqw (2D)
  w_present=0;
  q_present=0;
  if isa(s, 'iData') && ndims(s) == 2
    for index=1:2
      lab = lower(label(s,index));
      if any(strfind(lab, 'wavevector')) || any(strfind(lab, 'q')) || any(strfind(lab, 'Angs'))
        q_present=index;
      elseif any(strfind(lab, 'energy')) || any(strfind(lab, 'w')) || any(strfind(lab, 'meV'))
        w_present=index;
      end
    end
  end
  if ~w_present || ~q_present
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not seem to be an isotropic S(|q|,w) 2D object. Ignoring.' ]);
    return
  end

  % get the S(q,w) on a meshgrid (and possibly rebin/interpolate)
  sqw = meshgrid(s);
  
  % check if we need to transpose the S(q,w)
  if w_present==2 && q_present==1
    s = transpose(s);
  end
  
  % get axes
  w = sqw{1};
  q = sqw{2};

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
    Ki = max(max(q))/2;
    Ei = max(max(abs(w)));
    if SE2V*V2K*sqrt(Ei) < Ki % check which of q,w limits the range
      Ki = SE2V*V2K*sqrt(Ei);
    end
  end

  % compute Ei
  if isempty(Ki), return; end
  Vi = K2V*Ki;
  Ei = VS2E*Vi.^2;

  % handle a vector of incoming Ki
  if numel(Ki) > 1
    for index=1:numel(Ki)
      sqw = [ sqw feval(mfilename, s(index), Ei, angles) ];
    end
    return
  end

  % compute Kf,Ef,Vf,Vi,Ei,lambda for the dynamic range computation
  lambda = 2*pi./Ki;

  Ef = Ei-w;
  Vf = SE2V*sqrt(Ef);
  Kf = V2K*Vf;

  fprintf(1, '%s: incoming: Ei=%g [meV] Ki=%g [Angs-1] lambda=%g [Angs] Vi=%g [m/s]\n', ...
    mfilename, Ei, Ki, lambda, Vi);

  % find: cos theta = (ki2 + kf2 - q2)/(2ki kf) is NOT between -1 and 1. theta=diffusion angle
  costheta = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf);
  if isempty(angles)
    index= find(abs(costheta) > 1 | Ef <= 0);
  else
    cost = cosd(angles);
    index= find(costheta < min(cost(:)) | costheta > max(cost(:)) | Ef <= 0);
  end
  
  sqw(index)=0;
  sqw.DisplayName = strtrim([ sqw.DisplayName ' Ei=' num2str(Ei) ]);
  sqw = title(sqw, strtrim([ sqw.Title ' Ei=' num2str(Ei) ]));
  
