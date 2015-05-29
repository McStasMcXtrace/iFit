% waitfor(pid): wait for the end of a running process
%
% @Process/waitfor: wait or the end of a running Process. Its status
% output and error print-out are updated.
%
% input:
%   pid: a process object as created with Process('command')
%
% examples:
%   pid=Process('ping 127.0.0.1')
%   refresh(pid);
%   pid.stdout
%   waitfor(pid);
%
% See also: Process, Process/refresh, Process/exit
% Process is part of iFit http://ifit.mccode.org 

function pid=waitfor(pid)

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

pid.process.waitFor;
pid=refresh(pid);

if isempty(pid.terminationDate) || ~pid.terminationDate
  pid.terminationDate=now;
end
pid.isActive  = 0;
pid.exitValue = pid.process.exitValue;

if nargout == 0 && ~isempty(inputname(1))
  assignin('caller', inputname(1), pid);
end
