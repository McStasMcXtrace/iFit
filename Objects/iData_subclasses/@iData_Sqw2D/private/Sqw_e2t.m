function [sxt, schan] =Sqw_e2t(s, varargin)
% convert S(xx,w) to S(xx,t). Requires wavelength, chwidth, distance
%
% also returns the S(phi,channel)

  sxt = []; schan = []; 
  if isempty(s), return; end
  
  p = varargin2struct({'lambda' 'chwidth' 'distance' 'Ki' 'Ei'}, varargin, true);
  
  if isempty(p.lambda) && isfield(p, 'ei') && ~isempty(p.ei) && p.ei > 0
    p.lambda = sqrt(81.8042./p.ei);
  end
  if isempty(p.lambda) && isfield(p, 'ki') && ~isempty(p.ki) && p.ki > 0
    p.lambda = 2*pi./p.lambda;
  end
  
  if nargin < 2, lambda = []; end

  disp([ mfilename ': ' s.Tag ' ' s.Title ' Converting Axis 1 "' ...
    label(s,1) '": energy [meV] to time [sec].' ]);
  
  % search missing parameters
  [s,lambda,distance,chwidth,energy,wavevector] = Sqw_search_lambda(s);
  if isempty(p.lambda)    p.lambda = lambda; end
  if isempty(p.distance)  p.distance=distance; end
  if isempty(p.chwidth)   p.chwidth=chwidth; end
  
  % default missing parameters
  if isempty(p.lambda),   p.lambda   = 2.36; end
  if isempty(p.distance), p.distance = 4;    end

  SE2V = 437.393377;        % Convert sqrt(E)[meV] to v[m/s]
  V2K  = 1.58825361e-3;     % Convert v[m/s] to k[1/AA]
  K2V  = 1/V2K;
  VS2E = 5.22703725e-6;     % Convert (v[m/s])**2 to E[meV]
  Ki   = 2*pi/p.lambda;
  Vi   = K2V*Ki;
  Ei   = VS2E*Vi.^2;
  
  % compute final energy
  hw   = getaxis(s, 1);
  Ef   = Ei - hw;
  t_sample_detector = p.distance/Vi;      % this is the time for the elastic peak
  t    = Ei./Ef;
  t(t < 0) = nan;
  t    = t_sample_detector.*sqrt(t); % from sample
  
  if isempty(p.chwidth)
    t_valid = t(isfinite(t));
    p.chwidth = diff(unique(t_valid(:)));
    p.chwidth = mean(p.chwidth);
    disp([ mfilename ': ' s.Tag ' ' s.Title ' Using ChannelWidth=' num2str(p.chwidth) ' [s]' ]);
  end    % first channel
  
  dtdE    = t./(2.*Ef)*1e6;   % to be consistent with vnorm abs. calc. all times above calculated in sec.
  % kikf    = sqrt(Ei./Ef);   % The Kf/Ki correction is applied in Sqw2ddcs method
  sxt       = copyobj(s)./dtdE;
  
  sxt = iData(sxt); % make it a true iData
  setalias(sxt, 'time', t, 'Time of flight / sample [s]');
  setalias(sxt, 'ElasticPeakPosition', t_sample_detector, '[s] Elastic peak position');
  setalias(sxt, 'IncidentWavelength', p.lambda);
  setalias(sxt, 'ChannelWidth',       p.chwidth);
  setalias(sxt, 'Distance',           p.distance);
  setaxis(sxt, 1, 'time');
  
  sxt = commandhistory(sxt, 'e2t', s, p.lambda);
  sxt.Label = 'S(x, tof)';
  label(sxt, 0, [  'e2t' '(' label(s, 0) ')' ]);
  
  % generate a S(phi, channel) data set
  if nargout > 1
    if all(t(:) < 1) && p.chwidth > 1 % channel width in us
      p.chwidth = p.chwidth/1e6;
    end
    t = t / p.chwidth;
    % now we want EPP to be at t_sample_detector
    % t_elast = EPP*chwidth = t_sample_detector
    EPP       = t_sample_detector/p.chwidth;
    
    % make sure time channels start at 1
    deltaChan = min(t(:));
    t         = t-deltaChan+1;
    EPP       = EPP-deltaChan+1;
    
    schan = copyobj(sxt);
    setalias(schan, 'Channel', round(t), 'Time Channel [1]');
    setalias(schan, 'ElasticPeakPosition', EPP, '[chan] Channel for the elastic');
    setaxis(schan, 1, 'Channel');
    sxt.Label = 'S(x, tof_chan)';
    label(sxt, 0, [  'e2t[channel]' '(' label(s, 0) ')' ]);
  end
