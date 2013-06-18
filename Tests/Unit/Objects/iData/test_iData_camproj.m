function result=test_iData_camproj

  a=iData(peaks); 
  b=[ camproj(a) sum(a) prod(a) cumsum(a) cumprod(a) ];
  
  if length(b) == 5
    result =   'OK     iData camproj sum cumsum';
  else
    result =   'FAILED iData camproj sum cumsum';
  end
