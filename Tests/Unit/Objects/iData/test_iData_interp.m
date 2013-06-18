function result=test_iData_interp

  a=iData(peaks); b=interp(a, 'grid'); c=interp(a, 2);
  if all(size(a)*2 == size(c))
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
