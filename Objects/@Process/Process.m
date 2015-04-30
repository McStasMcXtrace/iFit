% pid=Process(command): spawn external processes
%
% The Process class allows to launch an external command asynchronously,
% inquire its output and state during execution, and terminate it.
%
% The command to launch can be given as a single string, or as a cellstr.
% When used with a single string for both the program and its arguments,
% note that the string is parsed using white space as the delimiter.
% This cellstr calling syntax allows to pass arguments containing space
% delimiters.
%
%   start an external command asynchronously
%     pid = Process(command)
%
%   determine if process is still running (active) and update
%   stdout/stderr 
%     refresh(pid)
%
%   wait for the end of the process (synchronous mode)
%     waitfor(pid)
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

function pid = Process(command)

if nargin == 0
  % create a new object
  pid.process = '';
  pid.command = '';
  pid.creationDate  = now;  % Creation date
  pid.terminationDate  = []; % end date
  pid.stdinStream   ='';
  pid.stdin = [];   % stores the stdout (yes!) from the process

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
