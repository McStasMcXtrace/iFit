function t_elast = Sqw_getEPP(s, t_present)
% Sqw_getEPP: compute the elastic peak position (EPP) and centre the time axis
  
  if nargin < 2, t_present=1; end
  
  t = getaxis(s, t_present);
  
  % in principle, the time structure is symmetric when integrated over wavevector
  % but not from time.
  EPP = Sqw_getT(s, {'ElasticPeakPosition' 'Elastic_peak_channel' 'Elastic_peak' 'Peak_channel' 'Elastic'});
  if ~isempty(EPP)
    t_elast0 = EPP;
  else t_elast0=[]; end
  
  % compute the EPP
  if ndims(s) > 1
    s_time = trapz(s, 2);       % compute time distribution, integrating angle
  else
    s_time = s;
  end
  [~,s_time_max]            = max(s_time, [],t_present);
  t_elast                   = t(s_time_max);  % EPP estimate from maximum
  dt = (max(t(:)) - min(t(:)))/10;  % select 10% around the maximum

  % get gaussian distribution around maximum
  s_time = xlim(s_time, t_elast+[-dt dt]);
  [s_time_std,s_time_centre]= std(s_time,t_present);
  if abs(t_elast - s_time_centre) < s_time_std && abs(t_elast - s_time_centre) < 1e-3
    t_elast = mean([t_elast s_time_centre]);  % improve accuray from gaussian std
  end
  disp([ mfilename ': Computed Elastic peak position (EPP) ' num2str(t_elast) ]);
  % check the EPP when stored/computed as channels
  if ~isempty(t_elast0) && t_elast0 && t_elast0 > 1 && t_elast > 1 ...
    if abs(t_elast0 - t_elast)/t_elast0 > 1e-2
      disp([ mfilename ': WARNING: ' s.Tag ' ' s.Title ' Elastic peak position mismatch. Using computed value.' ])
      disp([ '  The Elastic peak position from the file Parameters is ' num2str(t_elast0) ]);
    else t_elast = mean([t_elast0 t_elast]); end
  end
