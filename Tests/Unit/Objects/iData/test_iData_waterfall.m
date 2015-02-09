function result=test_iData_waterfall

  h= waterfall(iData(peaks));
  close(gcf);
  
  if ~isempty(h)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
