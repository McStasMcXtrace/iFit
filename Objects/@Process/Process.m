% pid=Process(command): spawn external processes
%
% The Process class allows to launch an external command in background,
% inquire its output and state during execution, and terminate it.
%
% The command to launch can be given as a single string, or as a cellstr.
% When used with a single string for both the program and its arguments,
% note that the string is parsed using white space as the delimiter.
% This cellstr calling syntax allows to pass arguments containing space
% delimiters.
%
%   start an external command asynchronously (in background)
%     pid = Process(command)
%
%   determine if process is still running (active) and update
%   stdout/stderr 
%     refresh(pid); status = pid.isActive;
%   or
%     status = isreal(pid);
%
%   wait for the end of the process (synchronous mode)
%     waitfor(pid)
%   or
%     while any(isreal(pid)); pause(1); display(pid); end
%
%   force process to end (kill)
%     exit(pid)
%
% input:
%   command: the external command line to launch, specified as a single 
%            string with any additional arguments, or as a cellstr
%            which first item is the commands, and subsequent items
%            are command parameters.
%            The command can also be given as an other Process object, which 
%            command is then re-used (to restart a similar process).
% output:
%   pid:     a Process object.
%
% examples:
%   pid=Process('ping 127.0.0.1')
%   pid=Process({'ping','127.0.0.1'})
%
% See also: Process/refresh, Process/exit
% Process is part of iFit http://ifit.mccode.org 
% (c) E.Farhi, ILL. License: EUPL.

function pid = Process(command)

if nargin == 0
  % create a new object
  pid.process = '';
  pid.command = '';
  pid.creationDate  = now;  % Creation date
  pid.terminationDate  = []; % end date
  pid.stdinStream   ='';
  pid.stdout = [];   % stores the stdout (yes!) from the process

  pid.stderrStream   = '';
  pid.stderr = [];   % stores the stderr from the process
  
  pid.exitValue=''; % only valid at end of process
  pid.isActive=0;
  pid.UserData=[];
  pid =  class(pid, 'Process');
else
  % convert input argument into object
  if ischar(command) % Process('command')
    % process: start: launch a command and initiate the process
    pid = Process;
    
    % See <https://docs.oracle.com/javase/6/docs/api/java/lang/Runtime.html>
    pid.process = java.lang.Runtime.getRuntime().exec(command);
    pid.command = command;
    pid.creationDate  = now;
    
    % I/O from process
    pid.stdinStream   = pid.process.getInputStream;
    pid.stderrStream  = pid.process.getErrorStream;

    % initial state and log gathering
    % pid = refresh(pid);
    
  elseif isa(command, 'Process')
    % process: restart a new one
    pid = Process(pid.command); % relaunch a similar process
  end
  
  % update object
  if isa(command, 'Process') && nargout == 0 && ~isempty(inputname(1))
    assignin('caller', inputname(1), pid);
  end
end

% ------------------------------------------------------------------------------
%creation:


%t = timer('TimerFcn',@TimerRefreshCmd, 
%  'StartFcn', { @TimerStartCmd, command },
%  'StopFcn', @TimerStopCmd,
%  'ExecutionMode', 'fixedSpacing',
%  Name, command,
%  'Period', 10.0);

%function TimerStartCmd(obj, event, command)

%  obj.UserData.process = java.lang.Runtime.getRuntime().exec(command);
%  obj.UserData.command = command;
%  obj.UserData.creationDate  = now;
%  obj.UserData.terminationDate=[];
%  obj.UserData.isActive  = 1;
%  obj.UserData.stdout = [];   % stores the stdout (yes!) from the process
%  obj.UserData.stderr = [];   % stores the stderr from the process
%  
%  % I/O from process
%  obj.UserData.stdinStream   = pid.process.getInputStream;
%  obj.UserData.stderrStream  = pid.process.getErrorStream;
%  
%  % update the Timer to hold the external process information
%  assignin('caller', inputname(1), obj);
%  
%function TimerStopCmd(obj, event)

%  obj.UserData.terminationDate=now;
%  obj.UserData.process.destroy;
%  pause(1) % wait a little for process to abort
%  obj.UserData.isActive  = 0;
%  obj.UserData.exitValue = obj.UserData.process.exitValue;
%  
%  % update the Timer to hold the external process information
%  assignin('caller', inputname(1), obj);
%  
%function TimerRefreshCmd(obj, event)

%  % first check if process is still valid. For this we try to get its exitValue
%  % which is only possible after completion and raises an error during execution.

%  try
%    obj.UserData.exitValue = obj.UserData.process.exitValue; % will raise error if process still runs
%    if isempty(obj.UserData.terminationDate) || ~obj.UserData.terminationDate
%      obj.UserData.terminationDate=now;
%    end
%    obj.UserData.isActive  = 0;
%    stop(obj);
%  catch
%    if isempty(obj.UserData.process) || isempty(obj.UserData.command)
%      obj.UserData.isActive  = 0;
%      stop(obj);
%    else
%      obj.UserData.isActive  = 1;
%    end
%  end

%  % then retrieve any stdout/stderr content (possibly in Buffer after end of process)
%  obj.UserData.stdout = strcat(obj.UserData.stdout, process_get_output(pid.stdinStream));
%  obj.UserData.stderr = strcat(obj.UserData.stderr, process_get_output(pid.stderrStream));
%  
%  % update the Timer to hold the external process information
%  assignin('caller', inputname(1), obj);
%  
%% private function to read a stream
%  function out = process_get_output(stream)
%  % process_get_output: get the stream (InputStream) content[private]
%  
%    out = '';
%    available = 0;
%    try
%      available = stream.available;
%    end
%    
%    % return when nothing to read or invalid stream
%    if available <= 0, return; end
%    
%    % Read the content of the stream.
%    %
%    % EPIC FAIL 1:
%    % The following method would be nice, but fails:
%    %   target_buffer = javaArray('java.lang.Byte', available);
%    %   stream.read(target_buffer, 0, available);
%    %
%    % Indeed, as matlab converts any Java array into a Matlab class, 
%    % the read method can not be used with additional parameters (signature not found)
%    % see http://www.mathworks.com/matlabcentral/answers/66227-syntax-for-call-to-java-library-function-with-byte-reference-parameter
%    %
%    % EPIC FAIL 2:
%    % using readLine from a BufferedReader stalls after a few iterations
%    % Reader       = java.io.InputStreamReader(stream);
%    % bufferReader = java.io.BufferedReader(Reader);
%    % readLine(bufferReader)
%    
%    % we use a for loop to read all bytes, one by one (sadly).
%    out = zeros(1, available);
%    for index=1:available
%      out(index) = read(stream);
%    end
%    
%    out = strrep(char(out), sprintf('\n\r'), sprintf('\n'));

%      
