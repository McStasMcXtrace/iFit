function result=test_iData_std

  a=iData(peaks);
  c=std(a); d=std(a,-1);
  
  if isreal(d) && ~isreal(c)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
