function result=test_iData_Sqw2D_deBosify

  s0=iData_Sqw2D('D2O_liq_290_coh.sqw.zip'); 
  T0 = 300;
  s =Bosify(symmetrize(s0), T0);
  
  s1 = deBosify(s);
  s1 = ylim(s1, [0 inf]);

  if sum(abs(s1 - s0),0) > 1e-3
    result = [ 'FAILED ' mfilename ];
  else
    result = [ 'OK     ' mfilename ];
  end
