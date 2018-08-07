function result=test_iData_Sqw2D_moments

  s0=iData_Sqw2D('D2O_liq_290_coh.sqw.zip'); 
  m0=moments(symmetrize(s0));
  
  s1=iData_Sqw2D('SQW_coh_lGe.nc'); s1.Temperature=1235; s1.weight=72.6;
  m1=moments(symmetrize(s1));
  
  r0=mean(m0(3)./m0(5));
  r1=mean(m1(3)./m1(5));

  % Wc and Wq are mostly proportional
  % 1st and 3rd moments are 0 for classical S(q,w)
  if abs(r0 - 4) > 1         || abs(r1 - 0.18) > .05 ...
  || abs(mean(m0(2))) > 1e-6 || abs(mean(m0(7))) > 1e-6 ...
  || abs(mean(m1(2))) > 1e-6 || abs(mean(m1(7))) > 1e-6
    result = [ 'FAILED ' mfilename ];
  else
    result = [ 'OK     ' mfilename ];
  end
