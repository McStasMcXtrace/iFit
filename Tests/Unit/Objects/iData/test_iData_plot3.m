function result=test_iData_plot3

  h= [ plot3(iData(peaks)) plot3(iData(flow)) ];
  
  if numel(h) == 2
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
