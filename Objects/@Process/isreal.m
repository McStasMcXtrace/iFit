% status = isreal(pid): read process status and output
%
% @Process/isreal: determines if the referred command is still running
% and collects the content of its printed output and error messages (stdout/stderr).
% The methods returns the state of the Process.
%
% input:
%   pid: a process object as created with Process('command')
% output:
%   status: true for active process, false when process is terminated.
%
% examples: waits for completion of process, equivalent to waitfor(Process(...))
%   pid=Process('ping 127.0.0.1')
%   while any(isreal(pid))
%     pause(1)
%   end
%     
%
% See also: Process/refresh, Process/exit
% Process is part of iFit http://ifit.mccode.org 

function status=isreal(pid)

% handle array of processes
if numel(pid) > 1
  status = zeros(size(pid))
  for index=1:numel(pid)
    status(index) = feval(mfilename, pid(index));
  end
else

  pid    = refresh(pid);
  status = pid.isActive;
end

% update the calling object
if ~isempty(inputname(1))
  assignin('caller', inputname(1), pid);
end

