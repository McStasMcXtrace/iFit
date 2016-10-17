function mifit_disp(message, log_only)
% [internal] mifit_disp: display message, and log it
  if nargin < 2, log_only=[]; end
  if isempty(log_only), log_only=false; end
  
  % handle input type: char, numeric, cellstr, ...
  if ~ischar(message) && isnumeric(message)
    message = cellstr(num2str(message));
  end
  
  if iscellstr(message)
    for index=1:numel(message)
      mifit_disp(message{index}, log_only)
    end
    return
  end
  
  % display
  if ~log_only
    if size(message,1) > 1
      disp(message);
    else
      disp([ 'mifit' ': ' message ]);
    end
  end
  
  % save to Log file 
  file = fullfile(prefdir, [ 'mifit' '.log' ]);
  fid = fopen(file, 'a+');
  if fid == -1, return; end
  fprintf(fid, '[%s] %s\n', datestr(now), message);
  fclose(fid);
