function result=test_miFit

  h = TextEdit;
  
  if isempty(h) || ~ishandle(h)
    result = [ 'FAILED ' mfilename ];
  else
    result = [ 'OK     ' mfilename ];
  end
  delete(h);



