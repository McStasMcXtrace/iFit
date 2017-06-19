function refresh_Process(pid)
  % check if a Process is still running. Collects its stdout/stderr.
  if length(pid) > 1
    % can not refresh an array.
    return
  end
  
  if ~isvalid(pid), return; end
  UserData = get(pid, 'UserData');
  try
    if isjava(UserData.process)
      UserData.exitValue = UserData.process.exitValue; % will raise error if process still runs
    end
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
    UserData.stdout = strcat(UserData.stdout, sprintf('\n'), toadd);
    if UserData.Monitor, disp(toadd); end
    
    toadd = process_get_output(UserData.stderrStream);
    UserData.stderr = strcat(UserData.stderr, sprintf('\n'), toadd);
    if UserData.Monitor, disp(toadd); end
    
  else
    UserData.stdinStream = [];
    UserData.stdoutStream= [];
  end
  % compute Duration
  UserData.Duration = etime(clock, datevec(UserData.creationDate));

  status = UserData.isActive;
  
  % when active, we execute the Callback
  if status
    Callback = UserData.TimerFcn;
    istop = exec_Callback(pid, Callback, 'refresh');
  end
  
  set(pid, 'UserData',UserData);
end
