function result=test_iData_contourf

  h= contourf(iData(peaks));
  close(gcf);
  if ~isempty(h)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
