function result=test_miFit

  % we must make sure the existing configuration is not lost
  % copy the .ini .mat and .log in tmp
  file = fullfile(prefdir, 'mifit.*');
  if ~isempty(dir(file))
    t = tempname;
    mkdir(t);
    copyfile(file, t);
    delete(file);
  else t = '';
  end
  
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
  if ~isempty(t)
    copyfile(fullfile(t, 'mifit.*'), prefdir);
    try
      rmdir(t, 's');
    end
  end
  
