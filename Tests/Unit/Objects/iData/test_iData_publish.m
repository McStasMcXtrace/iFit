function result=test_iData_publish

  a      = iData(peaks);
  tmp    = tempname; mkdir(tmp);
  result = [ 'OK     ' mfilename ]; failed = '';
  f      = save(a, fullfile(tmp,'test'), 'html'); 
  if isempty(f)
    result = [ 'FAILED ' mfilename ];
  end
  
  try; rmdir(tmp, 's'); end
  
  if ~isempty(failed)
    
  end
