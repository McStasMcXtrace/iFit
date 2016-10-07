function result=test_iData_event

  a=iData(rand(50,4));
  b=event(a);
  
  if ndims(a) == 2 && ndims(b) == 3
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
