function result=test_iData_caxis

  a=iData(peaks); 
  b=caxis(del2(a));
  
  result = [ 'OK     ' mfilename ];
