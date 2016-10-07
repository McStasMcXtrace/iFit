function result=test_iData_setaxis

  a=iData(peaks);
  setaxis(a, 1, 1:size(a,1) );
  b=rmaxis(a, 1);
  
  if ~isempty(getaxis(a, '1')) && isempty(getaxis(b, '1'))
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
