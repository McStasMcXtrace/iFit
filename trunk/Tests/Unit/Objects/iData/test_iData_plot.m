function result=test_iData_plot

  h= [ plot(iData(1:100)) plot(iData(peaks)) plot(iData(flow)) ];
  
  if numel(h) == 3
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
