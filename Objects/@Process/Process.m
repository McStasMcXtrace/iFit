classdef Process < timer
  % Process(command): starts a system command
    % 
    % pid = Process('command arguments ...')
    %
    % The Process class replaces the 'system' command. but is started asynchronously.
    % Matlab does not wait for the end of the Process to get back to interactive mode.
    % The stdout and stderr are collected periodically.  
    %
    % You can customize the Process with e.g. additional arguments such as:
    %   Process(..., 'TimeOut', value)  set a TimeOut (to kill Process after)
    %   Process(..., 'Period', value)   set the refresh rate in seconds (10 s).
    %   Process(..., 'Monitor', 0 or 1) flag to set the Process in silent/verbose mode
    %   Process(..., 'TimerFcn', @fcn)  execute periodically on refresh
    %   Process(..., 'StopFcn', @fcn)   execute when the Process is killed (stop/exit)
    %   Process(..., 'EndFcn', @fcn)    execute when the Process ends by itself
    %
    % The TimerFcn, StopFcn and EndFcn can be given as:
    %   * simple strings, such as 'disp(''OK'')'
    %   * a function handle with none to 2 arguments. The Callback will then 
    %     pass as 1st argument the Process object, and as 2nd the event
    %       in 'kill','timeout','end', or 'refresh'. 
    %     Example @(p,a)disp([ 'Process ' p.Name ': event ' a ])
    %   * the name of a function which takes none to 2 arguments. Same as above.
    % when a callback has a non-zero return value, it stops the Process.
    %
    % methods:
    %   refresh(pid)  force the pid to be refreshed, i.e check if it is running
    %                 and get its stdout/stderr.
    %   silent(pid)   set the process to silent mode (do not print stdout/stderr).
    %   verbose(pid)  set the process to verbose mode (do not print stdout/stderr).
    %   disp(pid)     display full Process information.
    %   pid           display short Process information. Same as display(pid).
    %   stdout(pid)   get the stdout stream from the Process (normal output).
    %   stderr(pid)   get the stderr stream from the Process (errors).
    %   exit(pid)     kill the Process (stop it). Same as stop(pid)
    %   delete(pid)   kill the Process and delete it from memory.
    %   waitfor(pid)  wait for the Process to end normally or on TimeOut.
    %
    % WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
    % This class derives from 'timer'. As such, it remains in memory after the 
    % Process has stopped. To remove it, use delete(pid), then clear(pid).
    % The UserData field of the object is used by the external command. In case
    % you need to store something, do NOT directly allocate obj.UserData, but 
    % rather store anything with:
    %    ud=get(pid,'UserData');
    %    ud.UserData = <your stuff>;
    %    set(pid, 'UserData', ud);
    % get your own data with:
    %    ud=get(pid,'UserData');
    %    <your stuff> = ud.UserData;
    %
    % Example:
    %   pid=Process('ping 127.0.0.1'); silent(pid);
    %   pause(5);
    %   exit(pid);
    %
    %   Copyright: Licensed under the BSD
    %              E. Farhi, ILL, France <farhi@ill.fr> Aug 2012, http://ifit.mccode.org
    
  properties
    % we store everything in the UserData of the timer object so that it fully
    % remains synchronized at all times, and is unique.
  end
  
  % --------------------------------------------------------------------------
  methods
    % the Process creator (initializer)
    function pid = Process(command, varargin)
      % Process(command): instantiate a process object and start the command
      
      % should add: Display = 1 => show stdout while it comes
      if nargin == 0, command = ''; end
      
      pid = pid@timer( ...
          'ExecutionMode', 'fixedSpacing', ...
          'Name', command, ...
          'UserData', [], ...
          'Period', 10.0);
          
      % Default properties
      UserData.process = [];       % Java RunTime object
      UserData.command = '';       % the command associated to the process
      UserData.creationDate = now;  % Creation date (start)
      UserData.terminationDate  = []; % end date
      UserData.stdinStream   ='';
      UserData.stdout=[];          % stores the stdout (yes!) from the process

      UserData.stderrStream='';
      UserData.stderr=[];          % stores the stderr from the process
    
      UserData.exitValue='';       % only valid at end of process
      UserData.isActive=0;
      UserData.TimerFcn='';        % executed everytime the refresh function is used
      UserData.StopFcn='';         % executed when process is stopped/killed
      UserData.EndFcn='';          % executed when process ends normally
      UserData.Monitor=1;
      UserData.UserData = [];      % if the User stores something it goes here with set/get overload
      UserData.TimeOut = [];
      UserData.Duration = 0;
      
      % search for specific options in varargin
      index = 1;
      while index <= numel(varargin)
        this = varargin{index};
        if index < numel(varargin), val = varargin{index+1}; 
        else val = []; end
        if ischar(this)
          switch lower(this)
          case 'monitor'
            if isempty(val), val = 1; end
            if ~isempty(val) && val, UserData.Monitor = val; end
          case 'silent'
            if isempty(val) , val = 0; end
            UserData.Monitor = ~val;
          case {'timerfcn'}
            if ~isempty(val), UserData.TimerFcn = val; end
          case {'endfcn','callback'}
            if ~isempty(val), UserData.EndFcn = val; end
          case {'StopFcn','killfcn'}
            if ~isempty(val), UserData.StopFcn = val; end
          case 'period'
            if ~isempty(val) && val > 0, set(pid, 'Period', val); end
          case 'userdata'
            if ~isempty(val), UserData.UserData = val; end
          case 'timeout'
            if ~isempty(val) && val>0, UserData.TimeOut = val; end
          end
        end
        index = index+1;
      end
          
      set(pid, 'TimerFcn',{@refresh_fcn, 'Refresh'}, ...
          'StopFcn', { @exit_fcn, 'Kill' });

      setUserData(pid, UserData);
          
      if ischar(command) && ~isempty(command)
        % create a timer object
        
        UserData.process = java.lang.Runtime.getRuntime().exec(command);
        UserData.command = command;
        UserData.creationDate  = now;
        
        % I/O from process
        UserData.stdinStream   = UserData.process.getInputStream;
        UserData.stderrStream  = UserData.process.getErrorStream;
        
        setUserData(pid,UserData);
        
        if UserData.Monitor
          disp([ datestr(now) ': Process ' get(pid,'Name') ' is starting.' ])
        end;
        start(pid);
      end
    end
    
    % --------------------------------------------------------------------------
    function refresh(pid)
      % Process/refresh(pid): poke a Process and update its stdout/stderr.
    
      refresh_Process(pid);

    end
    
    % --------------------------------------------------------------------------
    function ex = exit(pid)
      % Process/exit(pid): end/kill a running Process and/or return its exit value.
      % stop@timer(pid);
      if length(pid) == 1 && isvalid(pid) && any(strcmp(get(pid,'Running'),'on'))
        refresh_Process(pid);
      end
      stop(pid);
      ex = exit_Process(pid, 'kill');
      
    end
    
    % --------------------------------------------------------------------------
    function delete(pid)
      % Process/delete(pid): completely remove the Process from memory. 
      % The Process is killed. Its stdout/err/value are lost.
      
      if length(pid) ~= 1, return; end 
      if isvalid(pid) && strcmp(get(pid,'Running'),'on')
        exit(pid);
      end
      try
        delete@timer(pid); % remove the timer
      catch ME
        disp(ME.message)
      end

    end

    % i/o methods
    function s = stdout(pid)
      % Process/stdout(pid): return the standard output stream (stdout)
      if ~isvalid(pid), s=nan; return; end
      refresh(pid);
      UserData = getUserData(pid);
      s = UserData.stdout;
    end
    
    function s = stderr(pid)
      % Process/stderr(pid): return the standard error stream (stderr)
      if ~isvalid(pid), s=nan; return; end
      refresh(pid);
      UserData = getUserData(pid);
      s = UserData.stderr;
    end
    
    function s = isreal(pid)
      % Process/isreal(pid): return 1 when the process is running, 0 otherwise.
      if ~isvalid(pid), s=0; return; end
      refresh(pid);
      UserData = getUserData(pid);
      s = UserData.isActive;
    end
    
    function silent(pid)
      % Process/silent(pid): set the Process to silent mode.
      if ~isvalid(pid), return; end
      UserData = getUserData(pid);
      UserData.Monitor = 0;
      setUserData(pid, UserData);
    end
    
    function verbose(pid)
      % Process/verbose(pid): set the Process to verbose mode, which displays its stdout.
      if ~isvalid(pid), return; end
      UserData = getUserData(pid);
      UserData.Monitor = 1;
      setUserData(pid, UserData);
    end
    
    function waitfor(pid)
      % Process/waitfor(pid): wait for the Process to end normally or on TimeOut.
      % pressing Ctrl-C during the wait loop stops waiting, but does not kill the Process.
      if ~isvalid(pid), return; end
      period = get(pid, 'Period');
      while isreal(pid)
        pause(period);
      end
    end
    
    function stop(pid, action)
      stop@timer(pid);
      if nargin < 2, action='kill'; end
      exit_Process(pid,action);
    end
    
    function ud = getUserData(obj)
      ud = obj.jobject.UserData;
    end
    
    function setUserData(obj, ud)
      obj.jobject.UserData=ud;
    end

  end

end

% ------------------------------------------------------------------------------
% our timer functions
function refresh_fcn(obj, event, string_arg)

  refresh_Process(obj);
  UserData = get(obj, 'UserData');
  if ~UserData.isActive
    % Process has ended by itself or aborted externally
    disp([ mfilename ': Process ' get(obj,'Name') ' has ended.' ])
    stop(pid, 'end');
    
  elseif ~isempty(UserData.TimeOut) && UserData.TimeOut > 0 ...
    && etime(clock, datevec(UserData.creationDate)) > UserData.TimeOut
    stop(pid, 'timeout');
  end

end


function exit_fcn(obj, event, string_arg)
% called when the timer ends (stop)

  UserData = get(obj, 'UserData');
  if UserData.isActive && strcmp(get(obj,'Running'),'on')
    stop@timer(pid);
    exit_Process(obj, 'kill');  % kill
  end

end







