function [sxt, schan] =Sqw_e2t(s, lambda)
% convert S(xx,w) to S(xx,t). Requires wavelength, chwidth, distance
%
% also returns the S(phi,channel)

  sxt = []; schan = []; 
  if isempty(s), return; end
  if nargin < 2, lambda = []; end
  if isempty(lambda)
    [s,lambda,distance,chwidth] = Sqw_search_lambda(s);
  else
    [s,~,distance,chwidth] = Sqw_search_lambda(s);
  end

  disp([ mfilename ': ' s.Tag ' ' s.Title ' Converting Axis 1 "' ...
    label(s,1) '": energy [meV] to time [sec].' ]);
 
  if isempty(lambda),   lambda   = 2.36; end
  if isempty(distance), distance = 4; end

  SE2V = 437.393377;        % Convert sqrt(E)[meV] to v[m/s]
  V2K  = 1.58825361e-3;     % Convert v[m/s] to k[1/AA]
  K2V  = 1/V2K;
  VS2E = 5.22703725e-6;     % Convert (v[m/s])**2 to E[meV]
  Ki   = 2*pi/lambda;
  Vi   = K2V*Ki;
  Ei   = VS2E*Vi.^2;
  
  % compute final energy
  hw   = getaxis(s, 1);
  Ef   = Ei - hw;
  t_sample_detector = distance/Vi;      % this is the time for the elastic peak
  t    = Ei./Ef;
  t(t < 0) = nan;
  t    = t_sample_detector.*sqrt(t); % from sample
  
  if isempty(chwidth)
    t_valid = t(isfinite(t));
    chwidth = diff(unique(t_valid(:)));
    chwidth = mean(chwidth);
    disp([ mfilename ': ' s.Tag ' ' s.Title ' Using ChannelWidth=' num2str(chwidth) ' [s]' ]);
  end    % first channel
  
  dtdE    = t./(2.*Ef)*1e6;   % to be consistent with vnorm abs. calc. all times above calculated in sec.
  % kikf    = sqrt(Ei./Ef);   % The Kf/Ki correction is applied in Sqw2ddcs method
  sxt       = copyobj(s)./dtdE;
  
  sxt = iData(sxt); % make it a true iData
  setalias(sxt, 'time', t, 'Time of flight / sample [s]');
  setalias(sxt, 'ElasticPeakPosition', t_sample_detector, '[s] Elastic peak position');
  setalias(sxt, 'IncidentWavelength', lambda);
  setaxis(sxt, 1, 'time');
  
  sxt = commandhistory(sxt, 'e2t', s, lambda);
  sxt.Label = 'S(x, tof)';
  label(sxt, 0, [  'e2t' '(' label(s, 0) ')' ]);
  
  % generate a S(phi, channel) data set
  if nargout > 1
    if all(t(:) < 1) && chwidth > 1 % channel width in us
      chwidth = chwidth/1e6;
    end
    t = t / chwidth;
    % now we want EPP to be at t_sample_detector
    % t_elast = EPP*chwidth = t_sample_detector
    EPP       = t_sample_detector/chwidth;
    
    % make sure time channels start at 1
    deltaChan = min(t(:));
    t         = t-deltaChan+1;
    EPP       = EPP-deltaChan+1;
    
    schan = copyobj(sxt);
    setalias(schan, 'Channel', round(t), 'Time Channel [1]');
    setalias(schan, 'ElasticPeakPosition', EPP, '[chan] Channel for the elastic');
    setalias(schan, 'IncidentWavelength', lambda);
    setaxis(schan, 1, 'Channel');
    sxt.Label = 'S(x, tof_chan)';
    label(sxt, 0, [  'e2t[channel]' '(' label(s, 0) ')' ]);
  end
