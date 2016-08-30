function mifit_disp(message)
  % display message, and log it
  
  if size(message,1) > 1
    disp(message);
  else
    disp([ 'mifit' ': ' message ]);
  end
  file = fullfile(prefdir, [ 'mifit' '.log' ]);
  fid = fopen(file, 'a+');
  if fid == -1, return; end
  fprintf(fid, '[%s] %s\n', datestr(now), message);
  fclose(fid);
