function result = test_Treatment_Sqw

  Sqw   = iData_Sqw2D('SQW_coh_lGe.nc');
  Sqw2  = symmetrize(Sqw);
  SqwT  = Bosify(Sqw2, 1250);
  Sqw25 = dynamic_range(SqwT, 25);
  m=moments(SqwT, 72.64);
  XS    = scattering_cross_section(SqwT, 25); % 25 meV -> Ki=3.47
  file  = saveas(SqwT,'','mcstas');
  
  % total XS should be equal to \int(q S(q))/2Ki^2
  
  if abs(mean(m(1))-0.45) < 0.01 && abs(mean(m(2))-0.03) < 0.01 ...
  && abs(mean(m(3))-3.6) < 0.1    && abs(mean(m(4))-7.24) < 0.2 ...
  && abs(mean(m(5))-19.5) < 0.1  ...
  && abs(trapz(Sqw.q.*trapz(Sqw25)/2/3.47/3.47) - 0.47) < 0.01 ...
  && ~isempty(file)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end 
  delete(file);
