function ex=exit_Process(pid, action)
  % force to quit a running Process.
  if nargin < 2, action='end'; end

  if length(pid) > 1
    ex = nan;
    % can not kill an array.
    return
  end
  
  if ~isvalid(pid), ex=nan; return; end
  
  UserData = get(pid, 'UserData');
  if UserData.isActive
    % stop the timer but leaves the object. 
    if strcmp(get(pid,'Running'),'on'); stop(pid); end 
  else ex = UserData.exitValue; return;
  end
  
  UserData = get(pid, 'UserData');
  if isjava(UserData.process)
    UserData.process.destroy;
    pause(1) % wait a little for process to abort
  end
  
  if isjava(UserData.process)
    try
      UserData.exitValue = UserData.process.exitValue;
    catch
      % process is invalid (closed)
    end
  end
  ex = UserData.exitValue;
  if isempty(UserData.terminationDate) || ~UserData.terminationDate
    UserData.terminationDate=now;
  end
  UserData.process = [];
  UserData.isActive  = 0;
  % compute Duration
  UserData.Duration = etime(clock, datevec(UserData.creationDate));
  set(pid, 'UserData',UserData);
  
  refresh_Process(pid); % flush stdout/stderr

  % when active, we execute the Callback
  if strcmp(action,'kill') || strcmp(action,'timeout')
    Callback = UserData.StopFcn;
  elseif strcmp(action,'end')
    Callback = UserData.EndFcn;
  else Callback = '';
  end
  % display message
  if strcmp(action,'timeout')
    toadd = [ datestr(now) ': Process ' get(pid,'Name') ' has reached its TimeOut ' num2str(UserData.TimeOut) ' [s]' ];
    disp(toadd);
    UserData.stderr = strcat(UserData.stderr, sprintf('\n'), toadd);
    set(pid, 'UserData',UserData);
  elseif strcmp(action,'kill')
    toadd = [ datestr(now) ': Process ' get(pid,'Name') ' is requested to stop.' ];
    disp(toadd);
    UserData.stderr = strcat(UserData.stderr, sprintf('\n'), toadd);
    set(pid, 'UserData',UserData);
  end
  istop = exec_Callback(pid, Callback, action);
end
