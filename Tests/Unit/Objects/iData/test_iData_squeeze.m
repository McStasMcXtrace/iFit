function result=test_iData_squeeze

  a=iData(rand(5,1,3));
  c=squeeze(a);
  
  if length(size(a)) == 3 && length(size(c)) == 2
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
