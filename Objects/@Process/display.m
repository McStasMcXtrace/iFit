function d = display(s_in, name)
% d = display(s) : display Process object (from command line)
%
% @Process/display function to display Process object.
%   Used when no ';' sign follows a Process object in matlab.
% The return value may be catched as a string to display.  
%
% input:  s: object or array (Process) 
% output: d: string to display (char)
% ex:     'display(Process)' or 'Process'
%
% Version: $Revision$
% See also Process, Process/disp, Process/get

if nargin == 2 && ~isempty(name)
  iname = name;
elseif ~isempty(inputname(1))
  iname = inputname(1);
else
  iname = 'ans';
end

d = [ sprintf('%s = ',iname) ];

if numel(s_in) > 1
  d = [ d sprintf(' array [%s]',num2str(size(s_in))) ];
end
if isdeployed || ~usejava('jvm'), id='Process';
else           id='<a href="matlab:doc Process">Process</a> (<a href="matlab:methods Process">show methods</a>,<a href="matlab:doc(Process)">doc</a>)';
end
if length(s_in) == 0
    d = [ d sprintf(' %s: empty\n',id) ];
else
    if numel(s_in) == 1
      d = [ d sprintf(' %s:\n\n', id) ];
    else
      d = [ d sprintf(' %s:\n\n', id) ];
    end
    if numel(s_in) > 1
      d = [ d sprintf('Index ') ];
    end
    d = [ d sprintf('     [ID] [Command]                     [State] [output]\n') ];

    % now build the output string
    for index=1:numel(s_in)
      s = s_in(index);
      if length(s_in) > 1
        d = [ d sprintf('%5i ',index) ];                       % ID
      end
      c = char(s.process); if numel(c)>9, c=c((end-8):end); end
      d = [ d sprintf('%8s ', c) ];
      if iscellstr(s.command), c=sprintf('%s ', s.command{:});
      else c = char(s.command); end
      if numel(c)>30, c=[ c(1:27) '...' ]; end
      d = [ d sprintf('%30s ', s.command) ];                   % cmd;

      if s.isActive
        d = [ d 'Run    ' ];
      else
        d = [ d 'Stop   ' ];
      end
      d = [ d sprintf('%s\n', Process_display_out(s.stdin)) ];
    end
end

if nargout == 0
  fprintf(1,d);
end

% ------------------------------------------------------------------------------

function out = Process_display_out(str)
  if isempty(str), out=''; return; end
  lines = strread(str,'%s','delimiter','\n\r');
  out = sprintf('%s', lines{end});
  if numel(out) > 40, out = [ out(1:40) '...' ]; end
  
