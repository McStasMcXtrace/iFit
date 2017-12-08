function s=Sqw_t2e(s)
% convert S(xx,t) to S(xx,w). From lamp t2e and in5_t2e. Requires wavelength and/or distance

  if isempty(s), return; end
  [s,lambda,distance,chwidth] = Sqw_search_lambda(s);

  disp([ mfilename ': ' s.Tag ' ' s.Title ' Converting Axis 1 "' label(s,1) '": time [sec] to energy [meV].' ]);
  t = getaxis(s, 1);
  
  % check if the tof is given in channels
  if all(unique(diff(t(:))) == 1)
    % use ChannelWidth
    if ~isempty(chwidth) && chwidth
      t = t.*chwidth;
      disp([ mfilename ': the time axis is given in channels. Converting to time.' ]);
    else
      disp([ mfilename ': WARNING: ' s.Tag ' ' s.Title ' the time-of-flight Axis 1 is given in time channels.' ])
      disp('    This is probably NOT what you want. I will still try to use it as it is...')
      disp('    Define e.g. s.ChannelWidth=<channel width in time unit>');
    end
  end
  
  if  all(abs(t < .1))
    % probably in seconds (usually in the range 0 - 1e-2)
    % NOP: time is already OK
  elseif all(abs(t < 100))
    % probably in milli-seconds
    t = t/1000;
    disp('    Assuming time is in [ms].');
  elseif all(abs(t < 100000))
    % probably in micro-seconds
    t = t/1e6;
    disp('    Assuming time is in [us].');
  else
    disp([ mfilename ': WARNING: ' s.Tag ' ' s.Title ' the time-of-flight Axis 1 seems odd.' ])
    disp('    Check that the time-of-flight is defined as the time from the sample to the detector, in [s].')
  end
  
  % compute the elastic peak position (EPP) and centre the time axis
  % in principle, the time structure is symmetric when integrated over wavevector
  % but not from time.
  setaxis(s, 1, t);           % update time in [s]
  label(s, 1, 'time [s]');
  s_time = trapz(s, 2);       % compute time distribution, integrating angle
  [~,s_time_max]            = max(s_time, [],1);
  t_elast                    = t(s_time_max);  % EPP estimate from maximum
  dt = (max(t) - min(t))/10;  % select 10% around the maximum
  
  % get gaussian distribution around maximum
  s_time = xlim(s_time, t_elast+[-dt dt]);
  [s_time_std,s_time_centre]= std(s_time,1);
  if abs(t_elast - s_time_centre) < s_time_std && abs(t_elast - s_time_centre) < 1e-3
    t_elast = mean([t_elast s_time_centre]);  % improve accuray from gaussian std
  end
  
  % when angular axis, the estimate is wrong.
  % we should convert to [q,w], center energy range, and back to time.
  
  
  

  % search for the elastic peak position when lambda is not given
  if isempty(lambda) && ~isempty(t_elast) && ~isempty(distance)
    lambda   = t_elast./distance*3956.035;
  elseif ~isempty(lambda) && ~isempty(t_elast) && isempty(distance)
    distance = t_elast./lambda*3956.035;
  end
  
  % from there we need: t_elast, distance, lambda
  if isempty(lambda)   || lambda <= 0
    error([ mfilename ': could not determine incident neutron wavelength.' ]); end
  if isempty(t_elast)
    error([ mfilename ': could not determine elastic peak position.' ]); end
  if isempty(distance) || distance <= 0 
    error([ mfilename ': could not determine sample-detector distance.' ]); end
    
  % the time axis must be centered at t=0 on the sample
  SE2V = 437.393377;        % Convert sqrt(E)[meV] to v[m/s]
  V2K  = 1.58825361e-3;     % Convert v[m/s] to k[1/AA]
  K2V  = 1/V2K;
  VS2E = 5.22703725e-6;     % Convert (v[m/s])**2 to E[meV]
  
  Ki   = 2*pi/lambda;
  Vi   = K2V*Ki;
  Ei   = VS2E*Vi.^2;
  
  % center time to the sample
  t_sample_detector = distance/Vi;  % must correspond with EPP
  t                 = t -t_elast +t_sample_detector;  % shift time axis to EPP=sample-detector time
  
  Ef = Ei.*(t_sample_detector./t).^2;
  dtdE    = t./(2.*Ef)*1e6;   % to be consistent with vnorm abs. calc.
  kikf    = sqrt(Ei./Ef);     % all times above calculated in sec.
  hw = Ei - Ef;
  % Average energy scale:
  hw0 = mean(hw,2);

  dtdEkikf = dtdE.*kikf;
  s    = s.*dtdEkikf;
  setaxis(s, 1, hw0);
  label(s, 1, 'Energy transfer hw [meV]');
