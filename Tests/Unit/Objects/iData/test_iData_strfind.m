function result=test_iData_strfind

  a=iData('sv1850.scn')
  
  if ~isempty(strfind(a,'TITLE','case')) && ~isempty(strfind(a,'ILL TAS Data'))
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
