function result=test_iData_setalias

  a=iData(peaks);
  setalias(a,'NewStuff',20);
  setalias(a,'Signal', 'Data.Signal');
  setalias(a,'Signal', 'this.Data.Signal*2');
  b=rmalias(a, 'NewStuff');
  
  if isfield(a, 'NewStuff') && ~isfield(b, 'NewStuff');
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
