function [PID, command] = get_command(pid)
% poke for a PID/command in the current task list
%
% input:
%   pid: number or string to search in the processes
% output:
%   pid: the PID as a number, or empty when not found (int array)
%   command: the corresponding command (cellstr)

% pid     should be stored in Process.process property
% command should be stored in Process.command property

command = {}; PID = [];

if isobject(pid)
  ud = get(pid,'UserData');
  pid = ud.process;
  if isjava(pid)
    error([ mfilename ': This is/was a Java process.' ]);
  end
end

if isempty(pid), return; end
if isnumeric(pid) pid = pid(end); end

if ispc
  [response, tasks] = system('tasklist');
else
  [response, tasks] = system('ps -ef');
end

% now extract the PID and command
% split as lines
tasks = textscan(tasks, '%s', 'Delimiter',sprintf('\n')); tasks = tasks{1};

% should skip empty liens, those with 'PID' or '========'
if ispc
  % windows: processes are listed from line {4}. Lines {1:3} are comments.
  %   command is 1:25, PID is word following
  header= tasks{2};
  tasks = tasks(4:end);
else
  % Unix/MacOSX: processes are listed from line {2}. Lines {1} is comment.
  %   look for PID in task{1} (before). Look for CMD and get position (after)
  header  = tasks{1};
  tasks   = tasks(2:end);
end

% reduce the search list to PID/command matches
if isnumeric(pid), pid = num2str(pid); 
elseif ~ischar(pid)
  whos pid
  error([ mfilename ': PID specification should be a number or string' ]);
end
index = find(~cellfun(@isempty, strfind(tasks, pid)));
if isempty(index), command = {}; return; end
tasks = tasks(index);

% look for the PID and CMD in the header
index_PID = strfind(lower(header), ' pid ');
if isempty(index_PID)
  disp(header)
  warning([ mfilename ': Can not locate PID in process list.' ]);
else
  index_PID = index_PID(1);
end

for line = 1:numel(tasks)
  task = tasks{line};
  index_spc = find(task == ' ');
  if ~isempty(index_PID)
    % get the characters before the 'PID' token
    task_toPID = task(1:(index_PID+4)); % includes spaces around PID
    % search for the first space char before index_PID
    index_spc = index_spc(index_spc < index_PID); index_spc = index_spc(end);
    % get the PID
    PID = [ PID str2double(task(index_spc:(index_PID+4))) ];
  end
  % now get the 'command'
  command{end+1} = strtrim(task); % full line
end

