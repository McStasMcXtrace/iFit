function disp_Process(pid, name)
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
  else
    iname = 'ans';
  end
   
  if numel(pid) > 1
    eval([ iname ' = pid;' ])
    eval([ 'display(' iname ');' ]); % makes sure the variable name is sent to 'display'.
  else
    if ~isempty(inputname(1))
      refresh_Process(pid);
      assignin('caller', inputname(1), pid);
    end
    if isdeployed || ~usejava('jvm'), id='Process';
    else           id=[ '<a href="matlab:doc Process">Process</a> (<a href="matlab:methods Process">methods</a>,<a href="matlab:help Process">doc</a>,<a href="matlab:stdout(' iname ')">stdout</a>,<a href="matlab:exit(' iname ')">exit</a>)' ];
    end
    if ~isvalid(pid), return; end
    UserData = get(pid, 'UserData');
    if UserData.isActive, state='RUNNING'; else state='STOPPED'; end
    fprintf(1,'%s = %s object [%s]:\n',iname, id, state);

    s.Command      = UserData.command;
    s.creationDate = UserData.creationDate;
    s.terminationDate = UserData.terminationDate;
    s.exitValue    = UserData.exitValue;
    stdout = UserData.stdout;
    stderr = UserData.stderr;
    
    if isnumeric(s.creationDate),    s.creationDate=datestr(s.creationDate); end
    if isnumeric(s.terminationDate), s.terminationDate=datestr(s.terminationDate); end
    fprintf(1, '            process: %s\n', char(UserData.process));
    disp(s);
    % now display stdout/stderr tail
    if isdeployed || ~usejava('jvm') || ~usejava('desktop')
      fprintf(1, '             stdout: [%s char]\n', num2str(numel(stdout)));
    else
      fprintf(1, ['             <a href="matlab:stdout(' iname ')">stdout</a>: [%s char]\n'], num2str(numel(stdout)));
    end
    if numel(stdout), fprintf(1, Process_disp_out(stdout)); end
    if isdeployed || ~usejava('jvm') || ~usejava('desktop')
      fprintf(1, '             stderr: [%s char]\n', num2str(numel(stderr)));
    else
      fprintf(1,[ '             <a href="matlab:stderr(' iname ')">stderr</a>: [%s char]\n'], num2str(numel(stdout)));
    end
    if numel(stderr), fprintf(1, Process_disp_out(stderr)); end
  end
end

% ------------------------------------------------------------------------------
function out =  Process_disp_out(str)
  lines = strread(str,'%s','delimiter','\n\r');
  out = sprintf('                     %s\n', lines{1});
  if numel(lines) > 2, 
    out = [ out sprintf('                     ...\n') ];
  end
  if numel(lines) > 1, 
    out = [ out sprintf('                     %s\n', lines{end}) ];
  end
end
