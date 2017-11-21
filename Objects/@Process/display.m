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
  elseif length(s_in) >= 1
    % print header lines
    if length(s_in) == 1
      d = [ d sprintf(' %s:\n\n', id) ];
    else
      d = [ d id sprintf(' array [%s]',num2str(size(s_in))) sprintf('\n') ];
    end
    if length(s_in) > 1
      d = [ d sprintf('Index ') ];
    end
    d = [ d sprintf('     [ID] [Command]                     [State] [output]\n') ];

    % now build the output string
    for index=1:length(s_in)
      if length(s_in) > 1, d = [ d sprintf('%5i ', index) ]; end
      this = get_index(s_in, index);
      refresh_Process(this);
      UserData = get(this, 'UserData');
      if isjava(UserData.process)
        c = char(UserData.process);
      else c = num2str(UserData.process);
      end
      if numel(c)>9, c=c((end-8):end); end
      d = [ d sprintf('%8s ', c) ];
      if iscellstr(UserData.command), c=sprintf('%s ', UserData.command{:});
      else c = char(UserData.command); end
      if numel(c)>30, c=[ c(1:27) '...' ]; end
      d = [ d sprintf('%30s ', num2str(UserData.command)) ];                   % cmd;

      if UserData.isActive
        d = [ d 'Run    ' ];
      else
        d = [ d 'Stop   ' ];
      end
      if ~isempty(UserData.stderr), d=[ d '[ERR]' ]; end
      d = [ d sprintf('%s\n', Process_display_out(UserData.stdout)) ];
    end
  end

  if nargout == 0
    fprintf(1,d);
  end

end


% ------------------------------------------------------------------------------
function out = Process_display_out(str)
  if isempty(str), out=''; return; end
  lines = strread(str,'%s','delimiter','\n\r');
  if numel(lines) < 5
    out = sprintf('%s ', lines{:});
  else
    out = sprintf('%s ', lines{(end-4):end});
  end
  if numel(out) > 40, out = [ '...' out((end-35):end) ]; end
  out = deblank(out);
end

function obj = get_index(pid, index)
  S.type='()'; S.subs = {index};
  obj = subsref(pid, S);
end

