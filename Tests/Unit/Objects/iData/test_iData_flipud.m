function result=test_iData_flipud

  a=iData(peaks);
  b=flipud(a);
  c=fliplr(b);
  
  if abs(sum([a b c],0) - 836 < 1)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
