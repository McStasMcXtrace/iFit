function refresh_Process(pid)
  % check if a Process is still running. Collects its stdout/stderr.
  if length(pid) > 1
    % can not refresh an array.
    return
  end
  
  if ~isvalid(pid), return; end
  UserData = get(pid,'UserData');
  
  if isjava(UserData.process)
    pid = refresh_Process_java(pid);
  else
    pid = refresh_Process_external(pid);
  end
  UserData = get(pid,'UserData');
  
  % compute Duration
  UserData.Duration = etime(clock, datevec(UserData.creationDate));

  status = UserData.isActive;
  
  % when active, we execute the Callback
  if status
    Callback = UserData.TimerFcn;
    istop = exec_Callback(pid, Callback, 'refresh');
  end
  
  set(pid, 'UserData',UserData);


% ------------------------------------------------------------------------------
function pid = refresh_Process_java(pid)
  % test for a Java RunTime process (PID or command name)
  UserData = get(pid,'UserData');
  try
    UserData.exitValue = UserData.process.exitValue; % will raise error if process still runs
    if isempty(UserData.terminationDate) || ~UserData.terminationDate
      UserData.terminationDate=now;
    end
    UserData.isActive  = 0;
  catch ME
    % still running
    if isempty(UserData.process) || ~isjava(UserData.process)
      UserData.isActive  = 0;
    else
      UserData.isActive  = 1;
    end
  end
  
  % then retrieve any stdout/stderr content (possibly in Buffer after end of process)
  if isjava(UserData.stdinStream)
    toadd = process_get_output(UserData.stdinStream);
    UserData.stdout = sprintf('%s%s', UserData.stdout, toadd);
    if UserData.Monitor && ~isempty(toadd), disp([ get(pid,'Name') ': ' toadd ]); end
    
    toadd = process_get_output(UserData.stderrStream);
    UserData.stderr = sprintf('%s%s', UserData.stderr, toadd);
    if UserData.Monitor && ~isempty(toadd), disp([ get(pid,'Name') ': ' toadd ]); end
    
  else
    UserData.stdinStream = [];
    UserData.stdoutStream= [];
  end
  
  set(pid, 'UserData',UserData);

% ------------------------------------------------------------------------------
function pid = refresh_Process_external(pid)

  % test for an external process (PID or command name)
  
  [PID, command] = get_command(pid);
  UserData = get(pid,'UserData');
  
  if isempty(PID) && isempty(command)
    UserData.isActive  = 0;
    UserData.exitValue = 1;
  else
    UserData.isActive  = 1;
    UserData.exitValue = nan;
  end
  UserData.process = PID;
  if ~isempty(command)
    toadd = sprintf('%s\n', command{:});
    UserData.stdout = sprintf('%s%s', UserData.stdout, toadd);
    if UserData.Monitor && ~isempty(toadd), disp([ get(pid,'Name') ': ' toadd ]); end
  end
  set(pid, 'UserData',UserData);



