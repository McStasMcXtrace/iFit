function result = test_Treatment_Sqw

  Sqw   = load(iData, fullfile(ifitpath,'Data','SQW_coh_lGe.nc'));
  Sqw2  = Sqw_symmetrize(Sqw);
  SqwT  = Sqw_Bosify(Sqw2, 1250);
  Sqw25 = Sqw_dynamic_range(SqwT, 25);
  moments=Sqw_moments(SqwT, 72.64);
  XS    = Sqw_scatt_xs(SqwT, 25); % 25 meV -> Ki=3.47
  file  = Sqw_McStas(SqwT);
  
  % total XS should be equal to \int(q S(q))/2Ki^2
  
  if abs(mean(moments(1))-0.45) < 0.01 && abs(mean(moments(2))-1.26) < 0.01 ...
  && abs(mean(moments(3))-23.6) < 1    && abs(mean(moments(4))-5.4) < 0.1 ...
  && abs(mean(moments(5))-19.5) < 0.1  ...
  && abs(trapz(Sqw.q'.*trapz(Sqw25)/2/3.47/3.47) - 0.47) < 0.01 ...
  && ~isempty(file)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end 
  delete(file);
