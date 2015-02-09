function result=test_iData_scatter3

  h= scatter3(iData(peaks));
  close(gcf);
  
  if ~isempty(h)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
