function result=test_iData_end

  a = iData(peaks);
  
  b = zeros(a, 2,3);
  c = b(end);
  
  if all(c.Signal(:) == 1)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
