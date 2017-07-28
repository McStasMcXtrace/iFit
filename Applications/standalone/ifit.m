function ifit(varargin)
% Usage:  ifit [options] arguments commands 
%
%  iFit executable from command line.
%  Project page at <http://ifit.mccode.org>
%
%  arguments:
%    Any file name, including directories (imports all recursively).
%      Files can also be given as distant URLs (via http and ftp), and as compressed
%      content (zip, tar, gzip).
%      Script files (.m) are executed, other files are imported.
%    Any numerical value, including arrays given as e.g. '[ 1 2 ...]'.
%    The 'empty' argument is obtained from the '[]' argument
%    Any argument given between double quotes is used as a string argument.
%    Any argument given between simple quotes is evaluated as a an expression.
%
%  Once all arguments have been processed, they are stored in the 'this'
%    variable (cell array) for further use, and the program enters in interactive
%    mode except when the --exit argument is specified.
%  In interactive mode, any Matlab command can be entered if it spans on a single line.
%  In addition:
%    'verbose' and 'silent' commands toggle the echo of commands in the terminal.
%    'doc'     and 'help'   open the documentation
%    'run <script>'         executes a script, supporting multi-line
%    'propedit'             open the object browser GUI
%    'TextEdit'             open the script editor GUI, supporting multi-line
%    'mifit'                open the workspace GUI
%    
%  To execute multiple lines and control statements (if/while/for...), write them
%    in a script and type 'run <script>'.
%
%  options:
%  --exit or -e
%      exits - should be specified after all other commands/arguments
%  --save or -s or --save=FILE
%      save the workspace variables when commands have been executed.
%  --run=SCRIPT or -r=SCRIPT
%      executes the SCRIPT when starting.
%  -nodesktop
%      do not start the miFit GUI (prompt only) nor the TextEdit command window.
%  -desktop
%      start miFit and TextEdit command windows are opened
%
%  Examples:
%    ifit --save file1.*  subplot 

% Manual build:
% Better use the 'make.sh' script, or
% addpath(genpath('/home/farhi/svn/Matlab/iFit/trunk'))
% change all .m.org files into .m for the standalone
% mcc -m ifit -a /home/farhi/svn/Matlab/iFit/trunk
% buildmcr('.')

persistent ifit_options

if nargin == 1 && strcmp(varargin{1}, 'test')
  ifittest
  return
  
elseif ~isdeployed,
  disp([ mfilename ': you do not need the iFit Terminal. Just use Matlab as usual.' ])
  return; 
end
if ~isempty(ifit_options), return; end

inline_display_banner; % see inline below

ifit_options=inline_ifit_options(varargin{:});
this = {};

while ~exist('ifit_options') || ~isstruct(ifit_options) || ...
  (~strcmp(ifit_options.line, 'exit') && ~strcmp(ifit_options.line, 'return'))
  % restore internal ifit options if this was cleared
  if ~exist('ifit_options') || ~isstruct(ifit_options)
    ifit_options = inline_ifit_options;
  end
  ifit_options.line = strtrim(ifit_options.line);
  
  % handle specific commands (to override limitations from stand-alone)
  if strcmp(strtok(ifit_options.line,' '), 'doc')
    ifit_options.line = [ 'help' ifit_options.line(4:end) ];
    disp(ifit_options.line);
  end
  if strcmp(strtok(ifit_options.line,' ('),  'propedit')
    ifit_options.line = [ 'uiinspect' ifit_options.line(9:end) ];
  end
  if strcmp(strtok(ifit_options.line, ' ('), 'help')        % 'help' command ---------
    if length(ifit_options.line) > 4  % help <cmd>
      try
        ifit_options.line = inline_display_helpcommand(ifit_options.line); % returns empty line
      end
    else
      inline_display_help;            % see below (single 'help')
      ifit_options.line = 'doc(iData,''iFit''); disp('' '');';
    end
  elseif strncmp(ifit_options.line,'run ', 4) % 'run' command ------------------
    try
      ifit_options.line = inline_runscript(ifit_options.line);
    end
  elseif strncmp(ifit_options.line,'echo on', 7)  || strncmp(ifit_options.line,'verbose', 7)
    ifit_options.verbose = 1;
  elseif strncmp(ifit_options.line,'echo off', 8) || strncmp(ifit_options.line,'silent', 6)
    ifit_options.verbose = 0;
  end

  % argument is a file name/URL ?
  try
    mdir = dir(ifit_options.line);
  catch
    mdir = [];
  end
  if ~isempty(mdir) || ...
    any([ strncmp(ifit_options.line, {'file://','http://'},7) ...
          strncmp(ifit_options.line,  'ftp://', 6) ...
          strncmp(ifit_options.line,  'https://',8) ])
    if ~isempty(ifit_options.line)  % filenames have been dropped in the terminal
      if ~iscell(ifit_options.line), ifit_options.line = { ifit_options.line }; end
      ifit_options.line = ifit_options.line{1};
      this{end+1} = iData(ifit_options.line); % import them
      ifit_options.line = '';
      ans = this{end}
    end
  end
  
  % now do the work (evaluate what to do) --------------------------------------
  if ischar(ifit_options.line) && ~isempty(ifit_options.line), 
    if ifit_options.verbose
      disp(ifit_options.line)
    end
    try
      evalin('base',ifit_options.line); 
    catch ME
      disp(ME.message);
      if length(ifit_options.line) > 250
        ifit_options.line = [ ifit_options.line(1:250) ' ...' ];
      end
      inline_error(ifit_options.line);
      ifit_options.line = '';
    end
  end
  
  % collect next command to execute: from input arguments, or prompt
  if exist('varargin') == 1 && ~isempty(varargin) % from command line ----------
    % we clear the argument from the command line after reading it
    varargin = inline_cat_strings(varargin{:});
    ifit_options.line = varargin{1}; varargin(1) = [];

    % evaluate a char argument: can be a command, a "string", an 'expression'
    if isempty(ifit_options.line) || ~ischar(ifit_options.line), continue; end
    disp([ 'iFit:' num2str(ifit_options.index) '>> argument ' ifit_options.line ]);
    
    % specific case of imported arguments from the command line
    if (ifit_options.line(1)=='"' && ifit_options.line(end)=='"')
      % a "string" explicitly indicated as such
      ifit_options.line=ifit_options.line(2:(end-1));
      this{end+1} = ifit_options.line;
      ans = this{end}
      ifit_options.line = '';
    elseif (ifit_options.line(1)=='''' && ifit_options.line(end)=='''')
      % an 'expression' explicitly indicated as such
      ifit_options.line = ifit_options.line(2:(end-1));
      try
        [ifit_options.line] = evalin('base', ifit_options.line); % try to evaluate with a return argument
        this{end+1} = ifit_options.line;
        ans = this{end};
      catch
        try
          evalin('base', ifit_options.line) % evaluate without return argument
        catch ME
          disp('??? Error when evaluating expression argument:')
          disp(ifit_options.line)
          disp(ME.message)
          ifit_options.line = '';
        end
      end
      ifit_options.line = '';
    
    % some startup arguments known as commands
    elseif strcmp(ifit_options.line, '-nodesktop') || strcmp(ifit_options.line, '-nosplash')
      ifit_options.line = ''; % ignore these which are Matlab-desktop specific
      ifit_options.nodesktop=1;
    elseif strcmp(ifit_options.line, '-desktop')
      ifit_options.line = ''; % ignore these which are Matlab-desktop specific
      ifit_options.nodesktop=0;
      if feature('ShowFigureWindows')
        TextEdit({'% enter iFit/Matlab commands here', ...
          '% then use File/Evaluate', ...
          '% or the Run button'},'iFit commands');
        mifit;  % open mifit when no argin left and desktop/java on
    elseif strcmp(ifit_options.line, '--save') || strcmp(ifit_options.line, '-s')
      ifit_options.save='ifit.mat'; ifit_options.line = '';
    elseif strncmp(ifit_options.line, '--save=', 7)
      ifit_options.save=ifit_options.line(8:end); ifit_options.line = '';
    elseif strncmp(ifit_options.line, '--run=', 6)
      ifit_options.line=[ 'run ' ifit_options.line(7:end) ];
    elseif strncmp(ifit_options.line, '-r=', 3)
      ifit_options.line=[ 'run ' ifit_options.line(4:end) ];
    elseif strcmp(ifit_options.line, '-r')
      ifit_options.line=varargin{1}; varargin(1) = []; 
    elseif strcmp(ifit_options.line, '--exit') || strcmp(ifit_options.line, '-e')
      ifit_options.line='exit';
    elseif strcmp(ifit_options.line, '--help') || strcmp(ifit_options.line, '-h')
      inline_display_usage; % see below
      ifit_options.line = '';
      
    elseif strncmp(fliplr(ifit_options.line), fliplr('.m'), 2)
      % a script is given as argument : execute it
      ifit_options.line=[ 'run ' ifit_options.line ];
    elseif ~isempty(str2num(ifit_options.line))
      % numerical value(ifit_options.line) as a matrix
      this{end+1} = str2num(ifit_options.line);
      ans = this{end}; % we do not print the output, which may be BIG
      ifit_options.line = '';
    elseif ismethod(iData, ifit_options.line) || ismethod(iFunc, ifit_options.line) ...
            || (any(exist(ifit_options.line) == [ 2 3 5 6 ]) && isempty(dir(ifit_options.line)))
      % a known method/function (iData, iFunc, ...) but not a file name
      try
        ans = nargout(ifit_options.line);
        if nargout(ifit_options.line) > 1, 
          ans=cell(1,nargout(ifit_options.line));
          [ans{:}]          = builtin('feval',ifit_options.line, this{:});
          ifit_options.line = ans;
        else
          ifit_options.line = builtin('feval',ifit_options.line, this{:});
        end
      catch
        disp('??? Error when evaluating method:')
        disp(ifit_options.line)
        if ~isempty(this), disp(this); end
        disp(lasterr)
      end
      this{end+1} = ifit_options.line;
      ans = this{end}
      ifit_options.line = '';
    elseif strncmp(fliplr(ifit_options.line), fliplr('.desktop'), 8) || strncmp(fliplr(ifit_options.line), fliplr('.bat'), 4)
      % a desktop launcher is given as argument : read operator/command
      % from it and evaluate
      ifit_options.line = launcher_read(ifit_options.line);
      try
        this{end+1} = evalin('base', ifit_options.line); 
      catch
        this{end+1} = ifit_options.line;
      end
    else
      % argument is not an iFit/Matlab method, a command, a number, "string", 'expression', a script, a launcher
      % read data file and convert it to iData
      this{end+1} = iData(ifit_options.line);
      ifit_options.line = '';
      ans = this{end}
    end
    
  end

    ifit_options.index=ifit_options.index+1;
    
    if isempty(varargin) 
      if ~ifit_options.nodesktop && feature('ShowFigureWindows') && ifit_options.starting
        % open miFit and send imported objects there
        if numel(this) > 0
          inline_sendtomifit(this);
        end
      end
      if exist('this') && ~isempty(this) % last argument has just been processed
      
        % last file imported. Plot all imported data sets / functions
        if  isempty(ifit_options.line) % no command was given as last argument
          if iscell(this) && 0 < length(this(cellfun('isclass',this,'iData'))) && length(this(cellfun('isclass',this,'iData'))) <= 20
            figure('Name','iFit: imported data sets'); 
            subplot(this{cellfun('isclass',this,'iData')});
          elseif isa(this, 'iData') && numel(this) < 20
            figure('Name','iFit: imported data sets');
            subplot(this);
          end
          if iscell(this) && 0 < length(this(cellfun('isclass',this,'iFunc'))) && length(this(cellfun('isclass',this,'iFunc'))) <= 20
            figure('Name','iFit: imported models'); 
            subplot(this{cellfun('isclass',this,'iFunc')});
          elseif isa(this, 'iFunc') && numel(this) < 20
            figure('Name','iFit: imported models');
            subplot(this);
          end
        end
        
        if numel(this) > 1
          disp('Info: all imported arguments have been stored in cell array ''this''.');
          disp('      access them with e.g. this{1} ... this{end}');
          disp('      to get all models:    this(cellfun(''isclass'',this,''iFunc''))')
          disp('      to get all data sets: this(cellfun(''isclass'',this,''iData'')).')
        else
          this = this{1};
          disp('Info: imported argument has been stored in the array ''this''.');
        end
        disp('''this'' is:')
        disp(this);
        assignin('base', 'this', this);
        clear varargin
        this = {};
      end
      
    end
  else % not from varargin, from prompt ----------------------------------------
    ifit_options.line = input([ 'iFit:' num2str(ifit_options.index) ' ' ],'s');
    ifit_options.starting = 0;
  end
  if ~isempty(strtrim(ifit_options.line))
    ifit_options.index=ifit_options.index+1;
  end
end

% auto save ?
if ischar(ifit_options.save)
  save(ifit_options.save, '-v7.3'); % default format is Matlab/HDF (compressed)
end

close all
% close miFit and ResLibCal
try
  mifit File_Exit
  ResLibCal exit
  delete(findall(0,'Tag','TextEdit'))
end
disp([ '** Ending iFit on ' datestr(now) ])
disp(  '** Thanks for using iFit <ifit.mccode.org> **')
ifit_options = '';

% ------------------------------------------------------------------------------
%                             inline private functions
% ------------------------------------------------------------------------------

function inline_display_banner
  disp(' ');
  % banner from http://www.network-science.de/ascii/
  disp('                             _  ______ _  __ ');
  disp('                            (_)/ ____/(_)/ /_  (C) ILL');
  disp('                           / // /_   / // __/ ');
  disp('                          / // __/  / // /_   ');
  disp('                         /_//_/    /_/ \__/   ');
  disp(' ');
  disp('                          ** Welcome to iFit **');
  disp('                            <ifit.mccode.org>');
  disp('                E. Farhi, Institut Laue Langevin, France (c)');
  disp('                      Licensed under the EUPL V.1.1');
  disp(' ');
  disp(version(iData));
  disp(' ');
  disp([ '** Starting iFit on ' datestr(now) ])
  disp('   Type ''help'' to learn how to use this software.');
  disp('   iFit   help is fully available at <http://ifit.mccode.org>.')
  disp('   Matlab help is fully available at <http://www.mathworks.com/help/techdoc>.')
  disp('   Type ''exit'' or Ctrl-C to exit.');
  disp('** Applications (User Interfaces):');
  disp('   mifit (main GUI), rescal (neutron TAS), sqw_phonons (Phonons/DFT), ');
  disp('   sqw_spinw (SpinW, sw), TextEdit (editor/commands)');
  disp([ '** Example data files: ' fullfile(ifitpath ,'Data') ]);
  if ispc
    disp('WARNING: Windows: file names containing spaces, such as "My Documents" ')
    disp('         are not well supported. Rename files or move directories.')
  end
  if verLessThan('matlab', '7.11') 
    disp('WARNING: do NOT use accent characters at the prompt. This will CRASH. Use the TextEdit command to use a safe editor.');
  end
  disp(' ')

function inline_display_help
  disp(version(iData, 'contrib'))
  disp([ 'Built using Matlab ' version ])
  disp(' ');
  disp('Enter any Matlab/iFit command.');
  disp('      Use ''run script.m'' to execute a script from a file.');
  disp('      Control statements are allowed (for/while loops, tests, ');
  disp('        switch/case...) when they span on one line, or in scripts.');
  disp('Keys: Arrow-Up/Down  Navigate in command history.');
  disp('      Ctrl-C         Exit (same as ''exit'' or ''return'' commands)');
  disp('      Ctrl-u Ctrl-k  Delete to start/end of line');
  disp('      Ctrl-a Ctrl-e  Move cursor to start/end of line');
  disp('Help: Type ''doc(iData,''iFit'')'' to open local web pages');
  disp('      or see <ifit.mccode.org> and contact <ifit-users@mccode.org>.');
  disp('      Type ''help <ifit topic/method>'' to display specific help.')
  disp('In addition, you can use the commands:')
  disp('      verbose and silent    toggle the echo of commands in the terminal.')
  disp('      doc     and help      open the documentation')
  disp('      run <script.m>        executes a script, supporting multi-line')
  disp('      propedit              open the object browser GUI')
  disp('      TextEdit              open the script editor GUI, supporting multi-line')
  disp('      mifit                 open the workspace GUI')
  disp('      d=iData(''filename'') import a file')
  disp('      m=iFunc(''expr'');    create a model from an expression');
  disp('      fits(iFunc)           list of available models');
  disp('      fits(m,d)             fit a model to data');
  disp('Data and Models can be manipulated (+-/*...) using the Matlab syntax.');
  disp('Source code for this software is available at <ifit.mccode.org>.')
  disp('Matlab is a registered trademark of The Mathworks Inc.');
  disp('Matlab help is fully available at <http://www.mathworks.com/help/techdoc>.');
  
function inline_display_usage
  t = textwrap({ version(iData,'contrib') },80);
  fprintf(1, '%s\n', t{:});
  disp('Usage:  ifit [options] argruments commands ')
  disp(' ')
  disp('  iFit executable from command line.')
  disp('  Project page at <http://ifit.mccode.org>')
  disp(' ')
  disp('  arguments:')
  disp('    Any file name, including directories (imports all recursively).')
  disp('      Files can also be given as distant URLs (via http and ftp), and as compressed')
  disp('      content (zip, tar, gzip).')
  disp('      Script files (.m) are executed, other files are imported.')
  disp('    Any numerical value, including arrays given as e.g. ''[ 1 2 ...]''.')
  disp('    The ''empty'' argument is obtained from the ''[]'' argument')
  disp('    Any argument given between double quotes is used as a string argument.')
  disp('    Any argument given between simple quotes is evaluated as an expression.')
  disp(' ')
  disp('  Once all arguments have been processed, they are stored in the ''this''')
  disp('    variable (cell array) for further use, and the program enters in interactive')
  disp('    mode except when the --exit argument is specified.')
  disp('  In interactive mode, any Matlab command can be entered.')
  disp(' ')
  disp('  options:')
  disp('  --exit or -e')
  disp('      exits immediately after all execution of command line arguments.')
  disp('  --save or -s or --save=FILE')
  disp('      saves all variables when commands have been executed.')
  disp('  --run=SCRIPT or -r=SCRIPT')
  disp('      executes the SCRIPT when starting.')
  disp('  -nodesktop')
  disp('      does not start the miFit user interface')
  disp('  -desktop')
  disp('      make sure miFit and TextEdit command windows are opened')
  disp(' ')
  disp('  Examples:')
  disp('    ifit --save file1.*  subplot ')
  exit; % exit the application

function line = inline_display_helpcommand(line)
% handle 'help' and 'doc' commands.
  [t, line] = strtok(line); % remove first command 'help'
  line      = strtrim(line);
  if ~isempty(line)
    if isdeployed
      disp([ 'web(' line ')' ]);
      web(line);
    else
      help(line);
    end
    line = '';
  end  

function line = inline_runscript(line)
% load a script to execute
  line = strtrim(line(5:end)); % get script name
  if ~exist(line) && exist([ line '.m' ])
    line = [ line '.m' ];
  end
  if ~isempty(dir(line))
    [p,f,e]=fileparts(line);
    if isempty(p), p=pwd; end
    line = fullfile(p, [f e ]);
    disp([ 'Reading file ' line ]);
    line=fileread(line);
  else
    disp([ '??? Can not open script ' line ]);
  end
  
function inline_error(err)
  % display error trace
  disp('??? Error when evaluating expression:')
  disp(err)
  disp(lasterr)
  disp('Trace (first is where the error is detected):')
  for stack=getfield(lasterror,'stack')'
    if ischar(stack.name)
      fprintf('%s at line %i\n',stack.name,stack.line); 
    else
      disp(stack)
    end
  end

function varargin = inline_cat_strings(varargin)
  % test for arguments that start by " or ' and group them until we find similar ending char.
  if isempty(varargin), return; end
  v1 = varargin{1};
  if ~ischar(v1) || isempty(v1), return; end
  c = v1(1);
  if c == '"' || c == ''''
    for index=1:numel(varargin)
      % search for next argument which ends with " or '
      v2 = varargin{index};
      % concatenate to v1 if this is a string
      if index > 1 && ischar(v2)
        v1 = [ v1 ' ' v2 ];
        varargin{index} = ''; % this one has been moved into the 1st string 'v1'.
      end
      if ischar(v2) && v2(end) == c, break; end % we found the ending string token
    end
    % update string argument. Now contains {1:index}
    varargin{1} = v1;
  end
  
function inline_sendtomifit(this)
  % send 'this' to mifit
  h = mifit; % open miFit
  if isempty(this), return; end
  if iscell(this)
    for index=1:numel(this)
      mifit(this{index});
    end
  else
    mifit(this);
  end
  
function ifit_options=inline_ifit_options(varargin)
  ifit_options.verbose  =0;
  ifit_options.line     ='';     % the current line to execute
  ifit_options.index    =1;      % the index of the input
  this                  ={};     % the buffer from the command line
  ifit_options.save     =0;
  ifit_options.varargin =varargin;
  ifit_options.nodesktop=0;
  ifit_options.starting =1;
