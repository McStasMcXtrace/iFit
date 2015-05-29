% refresh(pid): read process status and output
%
% @Process/refresh: determines if the referred command is still running
% and collects the content of its printed output and error messages (stdout/stderr).
%
% input:
%   pid: a process object as created with Process('command')
%
% examples:
%   pid=Process('ping 127.0.0.1')
%   refresh(pid)
%   pid.stdout
%
% See also: Process/refresh, Process/exit
% Process is part of iFit http://ifit.mccode.org 

function pid=refresh(pid)

% handle array of processes
if numel(pid) > 1
  for index=1:numel(pid)
    pid(index) = feval(mfilename, pid(index));
  end
  if nargout == 0 && ~isempty(inputname(1))
    assignin('caller', inputname(1), pid);
  end
  return
end

% first check if process is still valid. For this we try to get its exitValue
% which is only possible after completion and raises an error during execution.

try
  pid.exitValue = pid.process.exitValue; % will raise error if process still runs
  if isempty(pid.terminationDate) || ~pid.terminationDate
    pid.terminationDate=now;
  end
  pid.isActive  = 0;
catch
  if isempty(pid.process) || isempty(pid.command)
    pid.isActive  = 0;
  else
    pid.isActive  = 1;
  end
end

% then retrieve any stdout/stderr content (possibly in Buffer after end of process)
pid.stdout = strcat(pid.stdout, process_get_output(pid.stdinStream));
pid.stderr = strcat(pid.stderr, process_get_output(pid.stderrStream));
  
if isa(pid, 'Process') && nargout == 0 && ~isempty(inputname(1))
  assignin('caller', inputname(1), pid);
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

      
