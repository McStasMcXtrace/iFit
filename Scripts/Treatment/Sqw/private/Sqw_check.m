function s = Sqw_check(s)
% Sqw_check: check if a 2D iData is aa S(q,w).
% input:
%   s: Sqw data set
%        e.g. 2D data set with w as 1st axis (rows), q as 2nd axis.
  
  % check if the data set is Sqw (2D)
  w_present=0;
  q_present=0;
  if isa(s, 'iData') && ndims(s) == 2
    for index=1:2
      lab = lower(label(s,index));
      if isempty(lab), lab=lower(getaxis(s, num2str(index))); end
      if any(strfind(lab, 'wavevector')) || any(strfind(lab, 'momentum')) || strcmp(strtok(lab), 'q') || any(strfind(lab, 'Angs'))
        q_present=index;
      elseif any(strfind(lab, 'energy')) || any(strfind(lab, 'frequency')) || strcmp(strtok(lab), 'w') || any(strfind(lab, 'meV'))
        w_present=index;
      end
    end
  end
  if ~w_present || ~q_present
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not seem to be an isotropic S(|q|,w) 2D object. Ignoring.' ]);
    s = [];
    return
  end

  % check if we need to transpose the S(q,w)
  if w_present==2 && q_present==1
    s = transpose(s);
  end
  
  % this is the weighting for valid data
  s(~isfinite(s)) = 0;
  
  % check 'classical'
  if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
    classical0 = s.classical;
    if numel(classical0) > 1
      s.classical = classical0(1);
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
    clear s_res s_opp
    
    % mean_log_ratio = mean(log_s_ratio,0);
    % std_log_ratio  = std(log_s_ratio,0);
    
    % compute the temperature from the Data
    % log_s_ratio should be a constant if S(q,w) contains Bose
    % then kT = w./log_s_ratio
    T         = w./log_s_ratio*11.6045; % 1 meV = 11.6045 K
    if any(isfinite(T))
      T         = mean(real(T),0);
    else T=NaN;
    end

    if isfinite(T) && T < -0.1
      % energy axis is reverted
      s{1} = -s{1};
      T    = -T;
    end

    % temperature stored ?
    T0        = Sqw_getT(s);
    if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
      classical0 = s.classical;
      if numel(classical0) > 1
        classical0 = classical0(1);
        s.classical = classical0;
      end
    else
      classical0 = [];
    end

    % log_s_ratio should be about 0 if S(q,w) is symmetric
    if isnan(T) || ~isfinite(T) || T <= 0.1 || T > 3000
      classical = 1;
      T         = Sqw_getT(s);
    else
      classical = 0;
    end

    % display warnings when mismatch is found
    if ~isempty(classical0) && ~isempty(classical) && classical0 ~= classical
      if   classical0, classical_str='classical/symmetric';
      else             classical_str='experimental/Bose/quantum'; end
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
      disp(['    indicates a ' classical_str ' S(|q|,w) 2D object, but the analysis of the data shows it is not.' ]);
    elseif isempty(classical0)
      setalias(s,'classical', classical);
    end

    if ~isempty(T0) && ~isempty(T) && ~(0.9 < T/T0 & T/T0 < 1.1)
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' S(|q|,w) 2D object from ' s.Source ]);
      disp(['    indicates a Temperature T=' num2str(T0) ' [K], but the analysis of the data provides T=' num2str(T) ' [K].' ]);
    end
    if isempty(T0) && ~isempty(T) && T > 0.1 && T < 3000
      disp([ mfilename ': INFO: Setting temperature T=' num2str(T) ' [K] for data set ' s.Tag ' ' s.Title ' S(|q|,w) 2D object from ' s.Source ]);
      s.Temperature = T;
    end
    
  end % energy axis has +/- 

  
