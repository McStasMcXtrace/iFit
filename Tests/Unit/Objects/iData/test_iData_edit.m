function result=test_iData_edit

  a = iData(peaks);
  h = edit(a);
  if ishandle(h)
    f = get(h, 'Parent');
  else h=0; f=[];
  end
  
  if ~isempty(h) && strcmp(get(h,'Type'),'uitable')
    result = [ 'OK     ' mfilename ];
  elseif ~isjava('jvm')
    result = [ 'OK     ' mfilename ' (skipped, no display)' ];
  else
    result = [ 'FAILED ' mfilename ];
  end
  
  delete(f);
