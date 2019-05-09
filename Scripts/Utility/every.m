function fixedDelay = every(t, command)
% EVERY execute a command every given period
%   T = EVERY(DELAY, COMMAND) executes the command after the given delay, in seconds.
%   The DELAY must be a single number, in seconds.
%   The COMMAND can be (see 'timer' documentation):
%     'code'               a single string (matlab code), 
%     @code(src,evnt, ...) a function handle with arguments source(timer) and 
%                          event type, and optional arguments.
%
%   To stop the execution, use STOP(T).
%
% Example: every(3, 'disp hello')
fixedDelay = [];
if nargin ~=2, return; end
if ~isnumeric(t) && ~isscalar(t)
  error([ mfilename ': invalid delay (' class(delay) '). Must be numeric/scalar.' ])
end
if (iscell(command) && numel(command) >= 1 && ...
  ~ischar(command{1}) && ~isa(command{1}, 'function_handle')) ...
  || (~ischar(command) && ~isa(command, 'function_handle') && ~iscell(command))
  error([ mfilename ': invalid command (' class(command) '). Must be char, cell or function_handle.' ])
end
if iscell(command) c=char(command{1}); else c =char(command); end
fixedDelay = timer('TimerFcn',command, 'ExecutionMode','fixedDelay', ...
    'Name', c, 'Period', t, 'StopFcn',@(s,e)delete(s));
start(fixedDelay);

% final message

disp([ mfilename ': starting command ' c ' with period ' num2str(t) ' seconds.' ])
