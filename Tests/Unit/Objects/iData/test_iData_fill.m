function result=test_iData_fill

  a=iData([ ifitpath 'Data/Monitor_GV*']); 
  
  b=hist(a); c=fill(b);
  
  if numel(find(isnan(b))) >= 0 && numel(find(isnan(c))) == 0
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
