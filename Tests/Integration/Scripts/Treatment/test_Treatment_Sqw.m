function result = test_Treatment_Sqw

  Sqw   = load(iData, fullfile(ifitpath,'Data','SQW_coh_lGe.nc'));
  Sqw2  = Sqw_symmetrize(Sqw);
  SqwT  = Sqw_Bosify(Sqw2, 1250);
  Sqw25 = Sqw_dynamic_range(SqwT, 25);
  moments=Sqw_moments(SqwT, 72.64);

  
  if abs(mean(moments(1))-0.45) < 0.01 && abs(mean(moments(2))-1.26) < 0.01 ...
  && abs(mean(moments(3))-23.6) < 1    && abs(mean(moments(4))-5.4) < 0.1 ...
  && abs(mean(moments(5))-19.5) < 0.1
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end 
