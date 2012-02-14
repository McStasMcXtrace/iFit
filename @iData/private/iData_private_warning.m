function iData_private_warning(a,b)
% function iData_private_warning
%
% iData_private_warning('enter');  stores the current warning level and unactivate 
%                                  further warning when level is higher than 1.
% iData_private_warning('exit');   restores the previous warning level state
% iData_private_warning('mfilename','MSG');

persistent warn;   % store the warning level

if nargin == 0
  a = 'enter';
end

if isempty(warn)  % create the persistent variable the first time
  warn.date   =clock;
  warn.level  =0;
  warn.lasterr='';
  lasterr('');
end

% if an error has occured more than 10 s before, we reset the warning state except for long operations
if ~strcmp(warn.lasterr, lasterr) && etime(clock, warn.date) > 10
  [st,i]=dbstack;
  long_exec={'fits','interp','ieval','cat','conv','fft','std','peaks','del2','gradient','jacobian','load'};
  long_flag=0;
  for index=1:length(long_exec)
    if any(strcmp(long_exec{index},{ st.name }))
      long_flag=1;
      break
    end
  end
  if ~long_flag
    warn.level  =0; % reset
  end
end
if ~strcmp(warn.lasterr, lasterr) && etime(clock, warn.date) < 1
  warn.level=warn.level-1;
end

% reset the warning level (e.g. when outdated/invalid)
if warn.level <=0 || strcmp(a, 'reset')
  warn.level = 1;
  warn.date  = clock;
  warn.caller= 'root';
  warn.lasterr='';
  warn.structs=[];
  warning('on','iData:setaxis');
  warning('on','iData:getaxis');
  warning('on','iData:get');
  warning('on','iData:subsref');
  lasterr('');
end

if strcmp(a, 'enter')
  % save the current warning state (level)
  % disable some warnings
  warn.structs = [ ...
    warning('off','iData:setaxis') warning('off','iData:getaxis') ...
    warning('off','iData:get')     warning('off','iData:subsref') ];
  warn.date = clock;
  warn.level=warn.level+1;
elseif strcmp(a, 'exit')
  % restore the previous warning state
  if ~isempty(warn.structs)
    warning(warn.structs);
  end
  warn.level=warn.level-1;
  warn.lasterr='';  % if we went there, all was fine: clean previous errors
  lasterr('');
  if warn.level > 50 % can not exceed too many levels
    warn.level  =0;  % reset (track of level goes wrong...)
  end
elseif nargin == 2
  warn.caller=a;
  b = [ 'iData/' a ': ' b ];  % MSG
  a = [ 'iData:' a ];         % ID
  warning(a,sprintf(b));
end

