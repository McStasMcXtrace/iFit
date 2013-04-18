function result = test_Math_10_interp
  a = iData(peaks(10))+2;
  b = interp(a,2);
  c = interp(a,1:.25:15,3:.25:12);
  a=iData([ ifitpath 'Data/Monitor_GV*']); 
  b=hist(a);
  result = 'OK  interp/hist';
