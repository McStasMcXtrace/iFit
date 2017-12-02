function s=Sqw_t2e(s)
% convert S(xx,t) to S(xx,w). From lamp t2e and in5_t2e. Requires wavelength and/or distance

  if isempty(s), return; end
  [s,lambda,distance,chwidth] = Sqw_search_lambda(s);

  disp([ mfilename ': ' s.Tag ' ' s.Title ' Converting Axis 1 "' label(s,1) '": time [sec] to energy [meV].' ]);
  t = s{1};
  
  % check if the tof is given in channels
  if all(unique(diff(t(:))) == 1)
    % use ChannelWidth
    if ~isempty(chwidth) && chwidth
      t = t.*chwidth;
    else
      disp([ mfilename ': WARNING: ' s.Tag ' ' s.Title ' the time-of-flight Axis 1 is given in time channels.' ])
      disp('    This is probably NOT what you want. I will still try to use it as it is...')
      disp('    Define e.g. s.ChannelWidth=<channel width in time unit>');
    end
  end
  
  if  all(t> 0 & t < .1)
    % probably in seconds (usually in the range 0 - 1e-2)
    % NOP: time is already OK
  elseif all(t > 0 & t < 100)
    % probably in milli-seconds
    t = t/1000;
    disp('    Assuming time is in [ms].');
  elseif all(t > 0 & t < 100000)
    % probably in micro-seconds
    t = t/1e6;
    disp('    Assuming time is in [us].');
  else
    disp([ mfilename ': WARNING: ' s.Tag ' ' s.Title ' the time-of-flight Axis 1 seems odd.' ])
    disp('    Check that the time-of-flight is defined as the time from the sample to the detector, in [s].')
    s = [];
    return;
  end
  
  % we compute the elastic peak position (EPP)
  telast = [];
  if isempty(lambda) || isempty(distance)
    s{2}   = t; % update time in [s]
    s_time = trapz(s, 2);
    [~,s_time_max]            = max(s_time, [],1);
    [s_time_std,s_time_centre]= std(s_time,1);
    if abs(t(s_time_max) - s_time_centre) < s_time_std && abs(t(s_time_max) - s_time_centre) < 1e-3
      telast = mean([t(s_time_max) s_time_centre]);
    else
      disp([ mfilename ': WARNING: ' s.Tag ' ' s.Title ' the time-of-flight Axis 1 elastic peak position seems odd.' ])
      disp('    Check that the time-of-flight is defined as the time from the sample to the detector, in [s].')
      disp('    Define e.g. s.Distance  =<sample-detector distance in m>');
      disp('    Define e.g. s.Wavelength=<lambda in Angs> or s.IncidentEnergy=<energy in meV>');
      s = [];
      return;
    end
  end

  % search for the elastic peak position when lambda is not given
  if isempty(lambda) && ~isempty(distance)
    lambda   = telast./distance*3956.035;
  elseif ~isempty(lambda) && isempty(distance)
    distance = telast./lambda*3956.035;
    if distance > 0
      setalias(s, 'Distance', distance, 'Sample-Detector distance [m]');
      disp([ mfilename ': ' s.Tag ' ' s.Title ' using <sample-detector distance> =' num2str(mean(distance(:))) ' [m]' ]);
    end
  else
    % we compute the elastic peak position (EPP)
    telast = distance .* lambda/3956.035;  % time from sample to detector elastic signal
  end
  
  if isempty(telast) || telast<= 0
    disp([ mfilename ': WARNING: ' s.Tag ' ' s.Title ' undefined sample-detector distance or wavelength.' ]);
    disp('    Define e.g. s.Distance=<sample-detector distance in m>');
    s = [];
    return;
  end
  
  Ei = 81.805./lambda^2;
  Ef = Ei.*(telast./t).^2;
  dtdE    = t./(2.*Ef)*1e6; % to be consistent with vnorm abs. calc.
  kikf    = sqrt(Ei./Ef);     % all times above calculated in sec.
  hw = Ei - Ef;
  % Average energy scale:
  hw0 = mean(hw,2);

  dtdEkikf = dtdE.*kikf;
  s    = s.*dtdEkikf;
  s{1} = hw0;
