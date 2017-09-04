function result=test_miFit

  % we must make sure the existing configuration is not lost
  % copy the .ini .mat and .log in tmp
  t = tempname;
  if isempty(dir(t))
    mkdir(t)
  end
  file = fullfile(prefdir, 'mifit.*');
  copyfile(file, t);
  delete(file);
  
  % now do the test
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

  % restore the mifit configuration files
  copyfile(fullfile(t, 'mifit.*'), prefdir);
  try
    rmdir(t, 's');
  end
