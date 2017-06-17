% exit(pid): kills a running process
%
% @Process/exit: requests the end of a running Process. Its status
% output and error print-out are updated.
%
% input:
%   pid: a process object as created with Process('command')
%
% examples:
%   pid=Process('ping 127.0.0.1')
%   refresh(pid);
%   pid.stdout
%   exit(pid);
%
% See also: Process/refresh, Process/exit
% Process is part of iFit http://ifit.mccode.org 

function pid=exit(pid)

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

pid=refresh(pid);
pid.process.destroy;
pause(1) % wait a little for process to abort

if isempty(pid.terminationDate) || ~pid.terminationDate
  pid.terminationDate=now;
  % execute DeleteFcn when defined
  if ~isempty(pid.DeleteFcn)
    if isa(pid.DeleteFcn, 'function_handle')
      n = nargin(pid.DeleteFcn);
      if n > 0
        feval(pid.DeleteFcn, pid);
      else
        feval(pid.DeleteFcn);
      end
    else
      eval(pid.DeleteFcn);
    end
  end
end

pid.isActive  = 0;
pid.exitValue = pid.process.exitValue;

if nargout == 0 && ~isempty(inputname(1))
  assignin('caller', inputname(1), pid);
end
