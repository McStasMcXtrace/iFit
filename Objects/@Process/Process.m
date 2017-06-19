classdef Process < timer
  properties
    % we store everything in the UserData of the timer object so that it fully
    % remains synchronized at all times, and is unique.
    UserData
  end
  
  % --------------------------------------------------------------------------
  methods
    % the Process creator (initializer)
    function pid = Process(command, varargin)
      % Process(command): starts a system command
      
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
      UserData.creationDate = [];  % Creation date (start)
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
          'StopFcn', { @exit_fcn, 'Killing process' });

      set(pid, 'UserData',UserData);
          
      if ischar(command) && ~isempty(command)
        % create a timer object
        
        UserData.process = java.lang.Runtime.getRuntime().exec(command);
        UserData.command = command;
        UserData.creationDate  = now;
        
        % I/O from process
        UserData.stdinStream   = UserData.process.getInputStream;
        UserData.stderrStream  = UserData.process.getErrorStream;
        
        set(pid, 'UserData',UserData);
        
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
      if strcmp(get(pid,'Running'),'on')
        refresh_Process(pid);
      end
  
      ex = exit_Process(pid, 'kill');
      
    end
    
    % --------------------------------------------------------------------------
    function delete(pid)
      % Process/delete(pid): completely remove the Process from memoty. 
      % The Process is killed. Its stdout/err/value are lost.
    
      exit(pid);
      try
        delete@timer(pid); % remove the timer
      catch ME
        disp(ME.message)
      end

    end

    % i/o methods
    function s = stdout(pid)
      % Process/stdout(pid): return the standard output stream (stdout)
      refresh(pid);
      UserData = get(pid, 'UserData');
      s = UserData.stdout;
    end
    
    function s = stderr(pid)
      % Process/stderr(pid): return the standard error stream (stderr)
      refresh(pid);
      UserData = get(pid, 'UserData');
      s = UserData.stderr;
    end
    
    function s = isreal(pid)
      refresh(pid);
      UserData = get(pid, 'UserData');
      s = UserData.isActive;
    end
    
    function silent(pid)
      UserData = get(pid, 'UserData');
      UserData.Monitor = 0;
      set(pid, 'UserData',UserData);
    end
    
    function verbose(pid)
      UserData = get(pid, 'UserData');
      UserData.Monitor = 1;
      set(pid, 'UserData',UserData);
    end
    
    function disp(pid)
      disp_Process(pid);
    end
    
    function display(pid)
      display_Process(pid);
    end
    
    
  end
  
  
    
end
      
% ------------------------------------------------------------------------------
% private function to read a stream
  function out = process_get_output(stream)
  % process_get_output: get the stream (InputStream) content[private]
  
    out = '';
    available = 0;
    try
      available = stream.available;
    end
    
    % return when nothing to read or invalid stream
    if available <= 0, return; end
    
    % Read the content of the stream.
    %
    % EPIC FAIL 1:
    % The following method would be nice, but fails:
    %   target_buffer = javaArray('java.lang.Byte', available);
    %   stream.read(target_buffer, 0, available);
    %
    % Indeed, as matlab converts any Java array into a Matlab class, 
    % the read method can not be used with additional parameters (signature not found)
    % see http://www.mathworks.com/matlabcentral/answers/66227-syntax-for-call-to-java-library-function-with-byte-reference-parameter
    %
    % EPIC FAIL 2:
    % using readLine from a BufferedReader stalls after a few iterations
    % Reader       = java.io.InputStreamReader(stream);
    % bufferReader = java.io.BufferedReader(Reader);
    % readLine(bufferReader)
    
    % we use a for loop to read all bytes, one by one (sadly).
    out = zeros(1, available);
    for index=1:available
      out(index) = read(stream);
    end
    
    out = strrep(char(out), sprintf('\n\r'), sprintf('\n'));

  end

% ------------------------------------------------------------------------------
% our timer functions
function refresh_fcn(obj, event, string_arg)

  refresh_Process(obj);
  UserData = get(obj, 'UserData');
  if ~UserData.isActive
    % Process has ended by itself or aborted externally
    disp([ mfilename ': Process ' get(obj,'Name') ' has ended.' ])
    exit_Process(obj, 'normal');
  elseif ~isempty(UserData.TimeOut) && UserData.TimeOut > 0 ...
    && etime(clock, datevec(UserData.creationDate)) > UserData.TimeOut
    disp([ mfilename ': Process ' get(obj,'Name') ' has reached its TimeOut ' num2str(UserData.TimeOut) ' [s]' ])
    exit_Process(obj, 'timeout');
  end

end


function exit_fcn(obj, event, string_arg)
% called when the timer ends (stop)

  UserData = get(obj, 'UserData');
  if UserData.isActive
    exit_Process(obj, 'kill');  % kill
  end

end

function refresh_Process(pid)
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

  status = UserData.isActive;
  
  % when active, we execute the Callback
  if status
    Callback = UserData.TimerFcn;
    istop = exec_Callback(pid, Callback, 'refresh');
  end
  
  set(pid, 'UserData',UserData);
end

function ex=exit_Process(pid, action)
  
  if nargin < 2, action='normal'; end
  
  UserData = get(pid, 'UserData');
  if UserData.isActive
    stop(pid); % stop the timer but leaves the object
  end
  
  UserData = get(pid, 'UserData');
  if isjava(UserData.process)
    UserData.process.destroy;
    pause(1) % wait a little for process to abort
  end
  
  UserData.isActive  = 0;
  if isjava(UserData.process)
    try
      UserData.exitValue = UserData.process.exitValue;
    catch
      % process is invalid (closed)
    end
  end
  UserData.process = [];
  set(pid, 'UserData',UserData);
  
  refresh_Process(pid);
  
  ex = UserData.exitValue;
  
  % when active, we execute the Callback
  if strcmp(action,'kill') || strcmp(action,'timeout')
    Callback = UserData.StopFcn;
  elseif strcmp(action,'normal')
    Callback = UserData.EndFcn;
  else Callback = '';
  end
  istop = exec_Callback(pid, Callback, action);
end

function istop = exec_Callback(pid, Callback, action)
  istop = 0; % failed ExternalFcns ignored
  if ~isempty(Callback) && (ischar(Callback) || isa(Callback, 'function_handle'))
    if isa(Callback, 'function_handle')
      nin = nargin(Callback);
      nout= nargout(Callback);
    else
      nin = 0; nout = 0;
    end
    
    if nin
      vars = { pid, action };
    end
    try
      if nout
        istop = feval(Callback, vars{1:nin});
      else
        eval(Callback);
      end
    catch ME
      disp(getReport(ME))
    end
    warning(action)
  end % Callback
end

function disp_Process(pid, name)
  % disp(s) : display Process object (details)
  %
  %   @Process/disp function to display Process object details
  %
  % input:  s: object or array (Process) 
  % ex:     'disp(Process)'
  %
  % Version: $Date$
  % See also Process, Process/refresh

  if nargin == 2 && ~isempty(name)
    iname = name;
  elseif ~isempty(inputname(1))
    iname = inputname(1);
  else
    iname = 'ans';
  end
   
  if numel(pid) > 1
    eval([ iname ' = pid;' ])
    eval([ 'display(' iname ');' ]); % makes sure the variable name is sent to 'display'.
  else
    if ~isempty(inputname(1))
      refresh_Process(pid);
      assignin('caller', inputname(1), pid);
    end
    if isdeployed || ~usejava('jvm'), id='Process';
    else           id='<a href="matlab:doc Process">Process</a> (<a href="matlab:methods Process">methods</a>,<a href="matlab:help Process">doc</a>)';
    end
    UserData = get(pid, 'UserData');
    if UserData.isActive, state='RUNNING'; else state='STOPPED'; end
    fprintf(1,'%s = %s object [%s]:\n',iname, id, state);

    s.Command      = UserData.command;
    s.creationDate = UserData.creationDate;
    s.terminationDate = UserData.terminationDate;
    stdout = UserData.stdout;
    stderr = UserData.stderr;
    
    if isnumeric(s.creationDate),    s.creationDate=datestr(s.creationDate); end
    if isnumeric(s.terminationDate), s.terminationDate=datestr(s.terminationDate); end
    fprintf(1, '            process: %s\n', char(UserData.process));
    disp(s);
    % now display stdout/stderr tail
    fprintf(1, '             stdout: [%s char]\n', num2str(numel(stdout)));
    if numel(stdout), fprintf(1, Process_disp_out(stdout)); end
    fprintf(1, '             stderr: [%s char]\n', num2str(numel(stderr)));
    if numel(stderr), fprintf(1, Process_disp_out(stderr)); end
  end
end

function d = display_Process(s_in, name)
  % d = display(s) : display Process object (from command line)
  %
  % @Process/display function to display Process object.
  %   Used when no ';' sign follows a Process object in matlab.
  % The return value may be catched as a string to display.  
  %
  % input:  s: object or array (Process) 
  % output: d: string to display (char)
  % ex:     'display(Process)' or 'Process'
  %
  % Version: $Date$
  % See also Process, Process/disp, Process/get

  if nargin == 2 && ~isempty(name)
    iname = name;
  elseif ~isempty(inputname(1))
    iname = inputname(1);
  else
    iname = 'ans';
  end

  d = [ sprintf('%s = ',iname) ];

  if numel(s_in) > 1
    d = [ d sprintf(' array [%s]',num2str(size(s_in))) ];
  end
  if isdeployed || ~usejava('jvm') || ~usejava('desktop'), id='Process';
  else           id=[ '<a href="matlab:doc Process">Process</a> (<a href="matlab:methods Process">methods</a>,<a href="matlab:help Process">doc</a>,<a href="matlab:disp(' iname ');">more...</a>)' ];
  end
  if length(s_in) == 0
      d = [ d sprintf(' %s: empty\n',id) ];
  else
      if numel(s_in) == 1 && ~isempty(inputname(1))
        refresh_Process(s_in);
        assignin('caller', inputname(1), s_in);
      end
      d = [ d sprintf(' %s:\n\n', id) ];
      if numel(s_in) > 1
        d = [ d sprintf('Index ') ];
      end
      d = [ d sprintf('     [ID] [Command]                     [State] [output]\n') ];

      % now build the output string
      for index=1:numel(s_in)
        s = s_in(index);
        UserData = get(s, 'UserData');
        if length(s_in) > 1
          d = [ d sprintf('%5i ',index) ];                       % ID
        end
        c = char(UserData.process); if numel(c)>9, c=c((end-8):end); end
        d = [ d sprintf('%8s ', c) ];
        if iscellstr(UserData.command), c=sprintf('%s ', UserData.command{:});
        else c = char(UserData.command); end
        if numel(c)>30, c=[ c(1:27) '...' ]; end
        d = [ d sprintf('%30s ', UserData.command) ];                   % cmd;

        if UserData.isActive
          d = [ d 'Run    ' ];
        else
          d = [ d 'Stop   ' ];
        end
        if ~isempty(UserData.stderr), d=[ d '[ERR]' ]; end
        d = [ d sprintf('%s\n', Process_display_out(UserData.stdout)) ];
      end
  end

  if nargout == 0
    fprintf(1,d);
  end

end
% ------------------------------------------------------------------------------
function out =  Process_disp_out(str)
  lines = strread(str,'%s','delimiter','\n\r');
  out = sprintf('                     %s\n', lines{1});
  if numel(lines) > 2, 
    out = [ out sprintf('                     ...\n') ];
  end
  if numel(lines) > 1, 
    out = [ out sprintf('                     %s\n', lines{end}) ];
  end
end

function out = Process_display_out(str)
  if isempty(str), out=''; return; end
  lines = strread(str,'%s','delimiter','\n\r');
  out = sprintf('%s', lines{end});
  if numel(out) > 40, out = [ out(1:40) '...' ]; end
  out = deblank(out);
end
