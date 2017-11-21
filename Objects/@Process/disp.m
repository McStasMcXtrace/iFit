function disp(s_in, name)
  % disp(s) : display Process object (details)
  %
  %   @Process/disp function to display Process object details
  %
  % input:  s: object or array (Process) 
  % ex:     'disp(Process)'
  %
  % Version: $Date$
  % See also Process, Process/refresh

  if nargin == 2 && ~isempty(name)
    iname = name;
  elseif ~isempty(inputname(1))
    iname = inputname(1);
  else
    iname = 'ans';
  end
   
  if length(s_in) > 1
    display(s_in, iname);
  else
    pid = s_in;
    refresh_Process(pid);
    if isdeployed || ~usejava('jvm'), id='Process';
    else           id=[ '<a href="matlab:doc Process">Process</a> (<a href="matlab:methods Process">methods</a>,<a href="matlab:help Process">doc</a>,<a href="matlab:stdout(' iname ')">stdout</a>,<a href="matlab:exit(' iname ')">exit</a>)' ];
    end

    if ~isvalid(pid), return; end
    UserData = get(pid, 'UserData');
    if UserData.isActive, state='RUNNING'; else state='STOPPED'; end
    fprintf(1,'%s = %s object [%s]:\n',iname, id, state);

    s.Command      = num2str(UserData.command);
    s.creationDate = UserData.creationDate;
    s.terminationDate = UserData.terminationDate;
    s.exitValue    = UserData.exitValue;
    stdout = UserData.stdout;
    stderr = UserData.stderr;
    
    if isnumeric(s.creationDate),    s.creationDate=datestr(s.creationDate); end
    if isnumeric(s.terminationDate), s.terminationDate=datestr(s.terminationDate); end
    if isjava(UserData.process)
      fprintf(1, '            process: %s\n', char(UserData.process));
    else
      fprintf(1, '            process: %s\n', num2str(UserData.process));
    end
    disp(s);
    % now display stdout/stderr tail
    if isdeployed || ~usejava('jvm') || ~usejava('desktop')
      fprintf(1, '             stdout: [%s char]\n', num2str(numel(stdout)));
    else
      fprintf(1, ['             <a href="matlab:stdout(' iname ')">stdout</a>: [%s char]\n'], num2str(numel(stdout)));
    end
    if numel(stdout), fprintf(1, '%s\n', Process_disp_out(stdout)); end
    if numel(stderr)
      if isdeployed || ~usejava('jvm') || ~usejava('desktop')
        fprintf(1, '             stderr: [%s char]\n', num2str(numel(stderr)));
      else
        fprintf(1,[ '             <a href="matlab:stderr(' iname ')">stderr</a>: [%s char]\n'], num2str(numel(stderr)));
      end
      fprintf(1, '%s\n', Process_disp_out(stderr)); 
    end
  end
end

% ------------------------------------------------------------------------------
function out =  Process_disp_out(str)
  if isempty(str), out=''; return; end
  lines = strread(str,'%s','delimiter','\n\r');
  if numel(lines) > 4
    n = numel(lines);
    out = sprintf('%s\n...\n%s\n%s\n', lines{1}, lines{n-1}, lines{n});
  else
    out = sprintf('%s\n', lines{:});
  end
  out = deblank(out);
end

function obj = get_index(pid, index)
  S.type='()'; S.subs = {index};
  obj = subsref(pid, S);
end
