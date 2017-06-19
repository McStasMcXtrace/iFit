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
  % Version: $Date$
  % See also Process, Process/disp, Process/get

  if nargin == 2 && ~isempty(name)
    iname = name;
  elseif ~isempty(inputname(1))
    iname = inputname(1);
  else
    iname = 'ans';
  end

  d = [ sprintf('%s = ',iname) ];

  if isdeployed || ~usejava('jvm') || ~usejava('desktop'), id='Process';
  else           id=[ '<a href="matlab:doc Process">Process</a> (<a href="matlab:methods Process">methods</a>,<a href="matlab:help Process">doc</a>,<a href="matlab:stdout(' iname ')">stdout</a>,<a href="matlab:exit(' iname ')">exit</a>,<a href="matlab:disp(' iname ');">more...</a>)' ];
  end
  
  if length(s_in) == 0
      d = [ d sprintf(' %s: empty\n',id) ];
  elseif length(s_in) == 1 && ~isvalid(s_in)
      d = [ d sprintf(' %s: invalid\n',id) ];
  elseif length(s_in) == 1
    if ~isempty(inputname(1))
      refresh_Process(s_in);
      assignin('caller', inputname(1), s_in);
    end
    d = [ d sprintf(' %s:\n\n', id) ];
    if length(s_in) > 1
      d = [ d sprintf('Index ') ];
    end
    d = [ d sprintf('     [ID] [Command]                     [State] [output]\n') ];

    % now build the output string

    UserData = get(s_in, 'UserData');
    if length(s_in) > 1
      d = [ d sprintf('%5i ',index) ];                       % ID
    end
    c = char(UserData.process); if numel(c)>9, c=c((end-8):end); end
    d = [ d sprintf('%8s ', c) ];
    if iscellstr(UserData.command), c=sprintf('%s ', UserData.command{:});
    else c = char(UserData.command); end
    if numel(c)>30, c=[ c(1:27) '...' ]; end
    d = [ d sprintf('%30s ', UserData.command) ];                   % cmd;

    if UserData.isActive
      d = [ d 'Run    ' ];
    else
      d = [ d 'Stop   ' ];
    end
    if ~isempty(UserData.stderr), d=[ d '[ERR]' ]; end
    d = [ d sprintf('%s\n', Process_display_out(UserData.stdout)) ];
  else
    d = [ d id sprintf(' array [%s]',num2str(size(s_in))) ];
    d = [ d sprintf('\n') ];
  end

  if nargout == 0
    fprintf(1,d);
  end

end


% ------------------------------------------------------------------------------
function out = Process_display_out(str)
  if isempty(str), out=''; return; end
  lines = strread(str,'%s','delimiter','\n\r');
  out = sprintf('%s', lines{end});
  if numel(out) > 40, out = [ out(1:40) '...' ]; end
  out = deblank(out);
end

