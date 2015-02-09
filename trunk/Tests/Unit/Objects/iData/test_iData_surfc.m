function result=test_iData_surfc

  h= surfc(iData(peaks));
  close(gcf);
  
  if ~isempty(h)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
