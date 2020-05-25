function result=test_iData_combine

  a=iData(peaks); a.Monitor=1;
  b=combine(a,a);
  
  if all(b.Monitor(:) == 2*a.Monitor(:)) && all(all(abs(sqrt(2)*a.Error - b.Error) < 1e-14)) ...
    && all(all(abs(2*a.Signal - b.Signal) < 1e-30))
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
  
