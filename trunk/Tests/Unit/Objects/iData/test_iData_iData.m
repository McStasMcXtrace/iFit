function result=test_iData_iData

  a=iData(peaks);
  
  if ~isempty(a) && isa(a, 'iData')
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
