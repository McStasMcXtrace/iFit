function result=test_iData_surfl

  h= surfl(iData(peaks));
  close(gcf);
  
  if ~isempty(h)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
