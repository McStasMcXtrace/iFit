function s = Sqw_check(s)
% Sqw_check: check if a 2D iData is a S(q,w).
%
% This routine can also convert automatically an input 
%         S(phi,t)  into S(phi,w)
%   and   S(phi,w)  into S(q,  w).
%
% conventions:
% omega = Ei-Ef = energy lost by the neutron
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
%
% input:
%   s: Sqw data set
%        e.g. 2D data set with w as 1st axis (rows), q as 2nd axis.
%
% Example: sqw=Sqw_check('SQW_coh_lGe.nc');
% (c) E.Farhi, ILL. License: EUPL.

  if nargin == 0, return; end
  
  if ~isa(s, 'iData'), s=iData(s); end
  
  % handle array of objects
  if numel(s) > 1
    sqw = [];
    for index=1:numel(s)
      sqw = [ sqw feval(mfilename, s(index)) ];
    end
    s(index)=iData; % free memory
    s = sqw;
    return
  end
  
  if isempty(s), return; end
  
  % check if the data set is Sqw (2D)
  w_present=0;
  q_present=0;
  a_present=0;
  t_present=0;
  alpha_present=0;
  beta_present=0;
  if isa(s, 'iData') && ndims(s) == 2
    for index=1:2
      lab = lower(label(s,index));
      def = getaxis(s, num2str(index));
      if ischar(def), lab = [ def ' ' lab ]; end
      if isempty(lab), lab=lower(getaxis(s, num2str(index))); end
      lab = strread(lab, '%s'); % split string into cell
      if strcmpm(lab, {'alpha','a'}) % strcmpm = multiple strcmpi is private below
        alpha_present=index;
      elseif strcmpm(lab, {'beta','b'})
        beta_present=index;
      elseif strcmpm(lab, {'wavevector','momentum','q','k','angs'})
        q_present=index;
      elseif strcmpm(lab, {'energy','frequency','w','e','mev'})
        w_present=index;
      elseif strcmpm(lab, {'time','sec','t','tof'})
        t_present=index;
      elseif strcmpm(lab, {'angle','deg','theta','phi'})
        a_present=index;
      end
    end
  end

  % conversions
  if alpha_present && beta_present && (~w_present || ~q_present)
    s = Sab_Sqw(s); % convert from S(alpha,beta) to S(q,w)
    return
  end
  
  % search for Sqw parameters for further conversions
  if ~isfield(s, 'parameters')
    s = Sqw_parameters(s, 'Sqw');
  end
  if ~w_present && t_present
    % convert from S(xx,t) to S(xx,w): t2e requires L2=Distance
    s = Sqw_t2e(s);
    s = Sqw_check(s);
    return
  end
  if ~q_present && a_present && w_present
    % convert from S(phi,w) to S(q,w)
    s = Sqw_phi2q(s);
    s = Sqw_check(s);
    return
  end
  if ~w_present || ~q_present
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
    disp('    does not seem to be an isotropic S(|q|,w) 2D object. Ignoring.');
disp(' ndims alpha  beta     q     w     t angle')
disp([ ndims(s) alpha_present beta_present q_present w_present t_present a_present ]);
    s = [];
    return
  end

  % check if we need to transpose the S(q,w)
  if w_present==2 && q_present==1
    s = transpose(s);
  end
  
  % this is the weighting for valid data
  s(~isfinite(s)) = 0;
  
  % check 'classical' and 'symmetric'
  if isfield(s,'classical') 
    classical0 = s.classical;
  else
    classical0 = [];
  end
  
  w  = s{1};
  % checks that we have 0<w<Ei for Stokes, and w<0 can be lower (anti-Stokes)
  if any(w(:) < 0) && any(w(:) > 0)
    % for experimental data we should get w1=max(w) < Ei and w2 can be as low as
    % measured (neutron gains energy from sample).
    
    
    w1 = max(w(:)); w2 = max(-w(:)); % should have w1 < w2
    if w1 > w2*2
      % we assume the measurement range is at least [-2*Ei:Ei]
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
      disp('    indicates that the energy range is mostly in the positive side.')
      disp('    Check that it corresponds with the neutron loss/sample gain Stokes side');
      disp('    and if not, revert energy axis with e.g. setaxis(s, 1, -s{1})');
    end
  end

  % can we guess if this is classical data ? get temperature ?
  w = s{1};

  if any(w(:) < 0) && any(w(:) > 0)
    % restrict the energy axis to the common +/- range
    w1 = max(w(:)); w2 = max(-w(:)); w_max = min([w1 w2]);

    if w1 ~= w_max || w2 ~= w_max
      s_res  = ylim(s, [-w_max w_max]); % restricted to [-w:w] range
    else
      s_res = s;
    end
    % get axes
    w = s_res{1};
    
    % we compare the s(q,w) and s(q,-w)
    s_opp = setaxis(s_res, 1, -w);
    s_opp = sum(s_opp,2); s_opp = sort(s_opp, 1);

    s_res = sum(s_res,2); s_res = sort(s_res, 1);

    % the ratio should be S(q,w)/S(q,-w) = exp(hw/kT)
    % so log(S(q,w)) - log(S(q,-w)) = hw/kT
    log_s_ratio = log(s_res) - log(s_opp);
    w = log_s_ratio{1};
    clear s_res s_opp
    
    % mean_log_ratio = mean(log_s_ratio,0);
    % std_log_ratio  = std(log_s_ratio,0);
    
    % compute the temperature from the Data
    % log_s_ratio should be a constant if S(q,w) contains Bose
    % then kT = w./log_s_ratio
    T         = w./log_s_ratio*11.6045; % 1 meV = 11.6045 K
    T         = T{0};
    if any(isfinite(T)) && (all(T(~isnan(T))>0.1) || all(T(~isnan(T))<0.1))
      T         = mean(real(T(~isnan(T))));
    else T=NaN;
    end

    if isfinite(T) && T < -0.1
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
      disp([ '    indicates a negative temperature T=' num2str(T) ' K. ' ]);
      disp(  '    Check the definition of the energy: Stokes=neutron losses energy=positive energy side');
      T = Sqw_getT(s);
    end

    % temperature stored ?
    T0        = Sqw_getT(s);

    % log_s_ratio should be about 0 if S(q,w) is symmetric
    classical = [];
    if ~isfinite(T) && ~isnan(T)
      classical = 1;
      T         = Sqw_getT(s);
    elseif T <= 0.1
      classical = 0;
    elseif T > 3000
      classical = 1;
    else 
    end

    % display warnings when mismatch is found
    if ~isempty(classical0) && ~isempty(classical) && classical0 ~= classical
      if   classical0, classical_str='classical/symmetric';
      else             classical_str='experimental/Bose/quantum/asymmetric'; end
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
      disp(['    indicates a ' classical_str ' S(|q|,w) 2D object, but the analysis of the data shows it is not.' ]);
    elseif isempty(classical0) && ~isempty(classical)
      setalias(s,'classical', classical);
    end

    if ~isempty(T0) && ~isempty(T) && ~isnan(T) && ~(0.9 < T/T0 & T/T0 < 1.1)
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' S(|q|,w) 2D object from ' s.Source ]);
      disp(['    indicates a Temperature T=' num2str(T0) ' [K], but the analysis of the data provides T=' num2str(T) ' [K].' ]);
    end
    if isempty(T0) && ~isempty(T) && ~isnan(T) && T > 0.1 && T < 3000
      disp([ mfilename ': INFO: Setting temperature T=' num2str(T) ' [K] for data set ' s.Tag ' ' s.Title ' S(|q|,w) 2D object from ' s.Source ]);
      s.Temperature = T;
    end
    
  end % energy axis has +/- 

% ------------------------------------------------------------------------------
function s = Sqw_phi2q(s)
% convert S(phi,w) to S(q,w). Requires wavelength

  [s,lambda] = Sqw_search_lambda(s);
  if isempty(s), return; end
  
  disp([ mfilename ': ' s.Tag ' ' s.Title ' Converting Axis 2: angle [deg] to wavevector [Angs-1].' ]);
  Ei    = 81.805/lambda^2;
  phi   = s{2}; % angle
  hw    = s{1};
  if isvector(hw) && isvector(phi)
    s = meshgrid(s);
    phi   = s{2}; % angle
    hw    = s{1};
  end
  % we use: cos(phi) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf);
  Ei = 81.805/lambda^2; Ki = 2*pi./lambda; 
  Ef = Ei - hw;         Kf = sqrt(Ef/2.0721);
  q  = sqrt(Ki.^2 + Kf.^2 - 2*cos(phi*pi/180).*Ki.*Kf);

  s = setalias(s, 'q', q, 'Wavevector [Angs-1]');
  s = setaxis(s, 2, 'q');
  
% ------------------------------------------------------------------------------
function s=Sqw_t2e(s)
% convert S(xx,t) to S(xx,w). From lamp t2e and in5_t2e. Requires wavelength and/or distance

  [s,lambda,distance,chwidth] = Sqw_search_lambda(s);
  if isempty(s), return; end

  disp([ mfilename ': ' s.Tag ' ' s.Title ' Converting Axis 1: time [sec] to energy [meV].' ]);
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
      disp('    Define e.g. s.Distance=<sample-detector distance in m>');
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
  else
    % we compute the elastic peak position (EPP)
    telast = distance .* lambda/3956.035;  % time from sample to detector elastic signal
  end
  
  if isempty(telast) || telast<= 0
    disp([ mfilename ': ' s.Tag ' ' s.Title ' undefined sample-detector distance or wavelength.' ]);
    disp('    Define e.g. s.Distance=<sample-detector distance in m>');
    s = [];
    return;
  end
  
  Ei = 81.805./lambda^2;
  Ef = Ei.*(telast./t).^2;
  dtdE    = t./(2.*Ef)*1e6; % to be consistent with vnorm abs. calc.
  kikf    = sqrt(Ei./Ef);            % all times above calculated in sec.
  hw = Ei - Ef;
  % Average energy scale:
  hw0 = mean(hw,2);

  dtdEkikf = dtdE.*kikf;
  s    = s.*dtdEkikf;
  s{1} = hw0;
  
% ------------------------------------------------------------------------------ 
function [s,lambda,distance,chwidth] = Sqw_search_lambda(s)
  % search for the wavelength etc in the object
  
  lambda = []; distance = []; chwidth=[];

  % no wavelength defined: search in object
  if isempty(lambda) && isfield(s, 'wavelength')
    lambda = mean(s.lambda);
    disp([ mfilename ': ' s.Tag ' ' s.Title ' using incident wavelength=' num2str(lambda) ' [Angs]' ]);
  end
  % search incident energy
  if isempty(lambda) && isfield(s, 'IncidentEnergy')
    energy = mean(s.IncidentEnergy);
    lambda = sqrt(81.805/energy);
    disp([ mfilename ': ' s.Tag ' ' s.Title ' using incident energy=' num2str(energy) ' [meV]' ]);
    setalias(s, 'IncidentEnergy', energy, 'Incident Energy [meV]');
  end
  % search incident wavevector
  if isempty(lambda) && isfield(s, 'IncidentWavevector')
    momentum = mean(s.IncidentWavevector);
    lambda=2*pi/momentum;
    disp([ mfilename ': ' s.Tag ' ' s.Title ' using incident wavevector=' num2str(momentum) ' [Angs-1]' ]);
  end
  if isempty(lambda)
    disp([ mfilename ': ' s.Tag ' ' s.Title ' undefined incident neutron wavelength/energy.' ]);
    disp('    Define e.g. s.Wavelength=<lambda in Angs> or s.IncidentEnergy=<energy in meV>');
  end
  
  if ~isfield(s, 'Wavelength')
    setalias(s, 'Wavelength', lambda, 'Incident Wavelength [Angs]');
  end
  
  % search for a sample-to-detector distance
  if isfield(s,'Distance'), 
    distance = s.Distance;
    disp([ mfilename ': ' s.Tag ' ' s.Title ' using <sample-detector distance> =' num2str(mean(distance(:))) ' [m]' ]);
  end
  
  % search for the Channel Width
  if isfield(s, 'ChannelWidth'), 
    chwidth = s.ChannelWidth;
    disp([ mfilename ': ' s.Tag ' ' s.Title ' using <channel width> =' num2str(mean(chwidth(:))) ' [time unit]']);
  end

% ------------------------------------------------------------------------------
function flag=strcmpm(str, words)
% multiple strcmp
%
% input:
%   str:   string or cellstr of tokens
%   words: string or cellstr of words to search

  flag = false;
  if ischar(str),   str=strread(str, '%s','delimiter',' ,; $()[]{}=|<>&"/\:"'''); end
  if ischar(words), words = strread(words, '%s'); end;
  for index=1:numel(words)
    if any(strcmpi(str, words{index}))
      flag = true; return;
    end
  end
