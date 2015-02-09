function result=test_iData_contour

  h= contour(iData(peaks));
  close(gcf);
  if ~isempty(h)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
