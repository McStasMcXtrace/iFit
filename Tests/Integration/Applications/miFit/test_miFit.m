function result=test_miFit

  h = mifit;
  result = [ 'FAILED ' mfilename ];
  if isempty(h) || ~ishandle(h)
    return
  end
  
  try
    mifit('Models_View_Parameters');
    mifit('File_Exit');
    result = [ 'OK     ' mfilename ];
  end


