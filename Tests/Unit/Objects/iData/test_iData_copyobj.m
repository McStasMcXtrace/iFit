function result=test_iData_copyobj

  b=copyobj(iData(peaks));
  
  result = [ 'FAILED ' mfilename ];
