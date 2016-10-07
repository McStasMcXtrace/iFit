function result=test_iData_circshift

  a=iData(peaks);
  shift = 5;
  b=circshift(a, 1, shift);
  [w,c]=std([a b], -1);
  
  if length(c) == 2 && abs(diff(c) - shift) < .1
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
