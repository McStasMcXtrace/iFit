function result=test_iData_edit

  a = iData(peaks);
  h = edit(a);
  f = get(h, 'Parent');
  
  if ~isempty(h) && strcmp(get(h,'Type'),'uitable')
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
  
  delete(f);
