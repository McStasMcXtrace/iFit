function ret=ifit(varargin)
% ifit [options] argruments: performs an iFit operation
%
% iFit executable from command line.
% Use ifit --help to get more help.

% create stand alone:
% remove Data directory
% mcc -m ifit -a /home/farhi/svn/Matlab/iFit
% unzip(buildmcr('.'))
% ./run_fit v76 ...

if nargin == 0
  varargin{1} = '--help';
end

% default behaviour
ret = 0;      % return code (all is fine)
% default options
options.property= '';   % initial options Property/filter (for single string search)
options.value   = '';   % initial value for the property=value filter
options.exit    = 0;    % exit directly ?
options.command ='plot';% default command
options.save    = 0;
options.verbose = 0;    % display all commands to output ?
do_command=-1;

% iterative loop, until all arguments have been processed
argv = {}; this=[];
for index=1:length(varargin)
  arg=varargin{index};
  if ischar(arg)  % should be mostly the case
    arg=strtrim(arg);
    if (arg(1)=='"' && arg(end)=='"')
      % a string explicitly indicated as such
      arg=arg(2:(end-1)); argv{end+1} = arg;
    elseif (arg(1)=='''' && arg(end)=='''')
      % an expression explicitly indicated as such
      arg=arg(2:(end-1)); argv{end+1} = arg;
      try
        arg=eval(arg);
      end
    elseif strcmp(arg, '--deploy')
      
      m = 'mcc -m ifit -a /home/farhi/svn/Matlab/iFit';
      m = [ m ' -a ' which('docopt') ];
      eval (m)
      return
    elseif (strncmp(arg, '-r=', 3) || strncmp(arg, '-r:', 3)) && length(arg) > 3
      % search for options
      options.command=arg(4:end);
      do_command=1;
    elseif (strncmp(arg, '--command=', 10) || strncmp(arg, '--command:', 10)) && length(arg) > 10
      options.command=arg(11:end); 
      do_command=1;
    elseif (strncmp(arg, '--filter=', 9) || strncmp(arg, '--filter:', 9)) && length(arg) > 9
      arg=arg(10:end);
      [property, value]=strtok(strtrim(arg),'=');
      if length(value > 1) value=strtrim(value(2:end)); else value=''; end
      % filter current objects
      options.property= property;
      options.values  = value;
    elseif strcmp(arg, '--exit')
      options.exit=1;
    elseif strncmp(arg, '--save=', 7)
      arg=arg(8:end);
      options.save=arg;
    elseif strcmp(arg, '--save')
      options.save=1;
    elseif strcmp(arg, '--verbose')
      options.verbose=1;
    elseif strcmp(arg, '--help') || strcmp(arg, '-h')
      t = textwrap({ version(iData,'contrib') },80);
      fprintf(1, '%s\n', t{:});
      disp(' ')
      NL = sprintf('\n');
      t= [ '  ifit [options] argruments: performs an iFit operation' NL ...
        ' ' NL ...
        '  iFit executable from command line.' NL ...
        '  Project page at <http://ifit.mccode.org' NL ...
        ' ' NL ...
        '  arguments:' NL ...
        '    Any file name, including directories (imports all recursively).' NL ...
        '      Files can also be given as distant URLs (via http and ftp), and as compressed' NL ...
        '      content (zip, tar, gzip).' NL ...
        '    Any numerical value, including arrays given as e.g. ''[ 1 2 ...]''.' NL ...
        '    The ''empty'' argument is obtained from the ''[]'' argument' NL ...
        '    Any argument given between double quotes is used as a string argument.' NL ...
        '    Any argument given between simple quotes is evaluated as a an expression.' NL ...
        ' ' NL ...
        '  options:' NL ...
        '  -r="eval string"' NL ...
        '  --command="eval string"' NL ...
        '      apply the specified command/operator on the imported objects, from the arguments.' NL ...
        '      The command specification can be any Matlab expression, including loops and tests.' NL ...
        '      Each command applies on objects imported before the --command is found.' NL ...
        '      More than one command can be given for sequential processing.' NL ...
        '      To call a script ''fun'', use ''--command=fun''' NL ...
        '      In commands, ''this'' is the array of imported data sets (iData array), and' NL ...
        '      ''argv'' refers to all defined arguments (as a cell)' NL ...
        '  --filter=keyword' NL ...
        '      searches in the imported objects for the keyword and selects these for the ' NL ...
        '      specified command (must preceed --command). Only applies when more than one ' NL ...
        '      object is processed.' NL ...
        '  --filter:property=value' NL ...
        '      value and select these for the specified command (must preceed --command).' NL ...
        '      Only applies when more than one object is processed.' NL ...
        '  --filter=OFF' NL ...
        '      removes the current filter and work consequently on all objects' NL ...
        '  --exit' NL ...
        '      do not wait for all displayed interfaces and figures to be closed by user ' NL ...
        '      before exiting, that is exits immediately after all execution of commands.' NL ...
        '  --save or --save=FILE' NL ...
        '      save all current objects when commands have been executed' NL ...
        '  --help or -h' NL ...
        '      display help' NL ...
        '  --doc' NL ...
        '      open documentation pages in a browser' NL ...
        ' ' NL ...
        '  Examples:' NL ...
        '  ifit file1.* "shifted" --save --command=plot ' ];
      disp(t);
      if ~isdeployed && exist('mcc')
        disp(' ');
        disp('  --deploy')
        disp('      deploy iFit as a standalone app')
      end
      return
    elseif strcmp(arg, '--doc')
      % open web page: docopt is needed, include it in the build
      p=doc(iData); d=docopt;
      if ispc,       
        system([ d ' ' p ]);
      elseif isunix, 
        system([ 'sh -l -c ' d ' ' p ]);
      end
    elseif ~isempty(str2num(arg))
      % numerical value(s) as a matrix
      arg=str2num(arg); argv{end+1} = arg;
    else
      % apparently not a file name ? -> store string as is
      % must not include filesep, not start with 'file','http','ftp',...
      % and must not exist in the path
      not_a_file =    isempty(strfind(arg, filesep)) ...
                   && isempty(strfind(arg, '\')) ...
                   && isempty(strfind(arg, '/')) ...
                   && isempty(strncmp(arg,{'file:','http:'},5)) ...
                   && isempty(strncmp(arg,'ftp:',4)) ...
                   && isempty(dir(arg));
      if ~not_a_file % probably a file name
        arg=iData(arg);
      end
      if ~isa(arg,'iData') || all(isempty(arg))
        arg=varargin{index}; % import did not work, store string directly
      end
      argv{end+1} = arg;
    end
    if do_command > 0
      % when no object has been defined, re-use the previous ones
      if ~any(cellfun('isclass', argv, 'iData')) && any(~isempty(this))
        argv{end+1} = this;
      end
      [ret, this] = ifit_do_command(options, argv{:});
      argv = {};
      do_command=0;
    end
  end
end % for

% when no command is specified, use default after parsing arguments
if do_command < 0 && ~isempty(argv)
  % no comand has been specified, use default
  [ret, this] = ifit_do_command(options, argv{:});
end

% save on exit ?
if numel(this)
  if isnumeric(options.save) && options.save
    save(this);
  elseif ischar(options.save)
    save(this, options.save);
  end
end

h = findobj('Type','figure');
% wait for all figures to be closed, except when --exit specified
if ~options.exit
  % now wait for all objects to be closed...
  for index=1:length(h)
    waitfor(h(index));
  end
end

% return: exit the command  
  
% ------------------------------------------------------------------------------
% private function to handle options until a command has been found
% then returns the previously defined options, and the iData array
function [ret, this]=ifit_do_command(options, varargin)

  % search for iData objects and gather them as a single array, in place
  % call this array 'this'
  this = []; first_iData=0; argv={};
  
  for index=1:length(varargin)
    if isa(varargin{index},'iData')
      if ~first_iData, first_iData=index; end
      that_one = varargin{index};
      if numel(that_one) > 1, that_one=reshape(that_one, 1, numel(that_one)); end
      this = [ this that_one ];
      argv{first_iData}=this;
    else
      argv{end+1} = varargin{index};
    end
  end
  
  % free memory
  clear index varargin that_one
  
  % filter objects 'this', if needed
  if ~isempty(options.property) && ~strcmpi(options.property, 'off') && numel(this) > 1
    to_keep=[];
    if isempty(options.value)
      for index=1:numel(this)
        if ~isempty(strfind(this(index), options.property))
          to_keep = [ to_keep this(index) ];
        end
      end
      if numel(to_keep), this=to_keep; end
    else
      this = findobj(this, options.property, options.value);
      if iscell(this), this=this{1}; end
    end
  end
  
  % display 'verbose'
  if options.verbose
    disp(options.command)
    for index=1:length(argv)
      this_arg=argv{index};
      if (isnumeric(this_arg) || ischar(this_arg)) && numel(this_arg) > 100
        this_arg = this_arg(1:100);
      end
      disp(this_arg);
    end
  end

  % evaluate the command
  try
    if isvarname(options.command)
      new_this=feval(options.command, argv{:});
      if isa(new_this,'iData'), this=new_this; end
    else
      eval(options.command);
    end
    ret=0;
  catch
    fprintf(1,'%s: Error while executing command %s\n', mfilename, options.command);
    ret=1;
  end
  
  if ~isa(this,'iData'), this=[]; end
  
