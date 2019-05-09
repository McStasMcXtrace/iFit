function singleShot = at(t, command)
% AT execute a command at given time
%   T = AT(DELAY, COMMAND) executes the command after the given delay, in seconds.
%   The DELAY can be a single number (in seconds), a date vector (6 
%   values from e.g. datevec and clock), a string such as '07-Apr-2008 23:00:00'.
%   The COMMAND can be (see 'timer' documentation):
%     'code'               a single string (matlab code), 
%     @code(src,evnt, ...) a function handle with arguments source(timer),
%                          event type, and optional arguments.
%
% Example: at(10, 'disp hello')
singleShot = [];
if nargin ~=2, return; end
if ischar(t), t = datenum(t)-clock; end
if isnumeric(t) && numel(t) == 6
  t = datenum(t)-now;
end
if ~isnumeric(t) && ~isscalar(t)
  error([ mfilename ': invalid delay (' class(delay) '). Must be char(date) or numeric/scalar.' ])
end
if (iscell(command) && numel(command) >= 1 && ...
  ~ischar(command{1}) && ~isa(command{1}, 'function_handle')) ...
  || (~ischar(command) && ~isa(command, 'function_handle') && ~iscell(command))
  error([ mfilename ': invalid command (' class(command) '). Must be char, cell or function_handle.' ])
end
if iscell(command) c=char(command{1}); else c =char(command); end
singleShot = timer('TimerFcn',command, 'ExecutionMode','singleShot', ...
    'Name', c, 'StartDelay', t, 'StopFcn',@(s,e)delete(s));
start(singleShot);

% final message
if iscell(command) c=char(command{1}); else c =char(command); end
disp([ mfilename ': starting command ' c ' in ' num2str(t) ' seconds.' ])
