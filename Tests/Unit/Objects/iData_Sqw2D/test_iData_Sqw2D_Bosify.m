function result=test_iData_Sqw2D_Bosify

  s0=iData_Sqw2D('SQW_coh_lGe.nc'); 
  T0 = 1235;
  s =Bosify(symmetrize(s0), T0);
  
  % now compute the temperature from the Bose factor (code from Sqw_getT)
  w = getaxis(s,1);
  % restrict the energy axis to the common +/- range
  w1 = max(w(:)); w2 = max(-w(:)); w_max = min([w1 w2]);

  if w1 ~= w_max || w2 ~= w_max
    s_res  = ylim(s, [-w_max w_max]); % restricted to [-w:w] range
  else
    s_res = s;
  end
  % get axes
  w = getaxis(s_res,1);
  
  % we compare the s(q,w) and s(q,-w)
  s_opp = setaxis(s_res, 1, -w);
  s_opp = sum(s_opp,2); s_opp = sort(s_opp, 1);

  s_res = sum(s_res,2); s_res = sort(s_res, 1);

  % the ratio should be S(q,w)/S(q,-w) = exp(hw/kT)
  % so log(S(q,w)) - log(S(q,-w)) = hw/kT
  wrn = warning;
  warning off;
  log_s_ratio = log(s_res) - log(s_opp);
  warning(wrn);
  w = getaxis(log_s_ratio,1);
  clear s_res s_opp
  
  % mean_log_ratio = mean(log_s_ratio,0);
  % std_log_ratio  = std(log_s_ratio,0);
  
  % compute the temperature from the Data
  % log_s_ratio should be a constant if S(q,w) contains Bose
  % then kT = w./log_s_ratio
  T         = w./log_s_ratio*11.6045; % 1 meV = 11.6045 K
  T         = getaxis(T,0);
  if any(isfinite(T)) && (all(T(~isnan(T))>0.1) || all(T(~isnan(T))<0.1))
    T         = mean(real(T(~isnan(T))));
  else T=0;
  end 
  
  if abs(T - T0) > 10
    result = [ 'FAILED ' mfilename ];
  else
    result = [ 'OK     ' mfilename ];
  end
