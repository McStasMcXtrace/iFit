function ifit_console(varargin)
% mcc -m ifit_console -a /home/farhi/svn/Matlab/iFit
% buildmcr('.')

disp(' ');
% banner from http://www.network-science.de/ascii/
disp('                      ____  ,                 _______ __ ');
disp('                     /  _/ // __    _ _      / ____(_) /_');
disp('                     / /     / /   / /      / /_  / / __/');
disp('                   _/ /     / /___/ /___   / __/ / / /_  ');
disp('                  /___/    /_____/_____/  /_/   /_/\__/  ');
disp(' ');
disp('                          ** Welcome to iFit **');
disp('                            <ifit.mccode.org>');
disp('                E. Farhi, Institut Laue Langevin, France (c)');
disp(' ');
disp([ '** Starting iFit on ' datestr(now) ])
disp(' ')

ifit_options.line     ='help'; % the current line to execute
ifit_options.index    =1;      % the index of the input
this                  ={};     % the buffer from the command line
ifit_options.save     =0;
ifit_options.exit     =0;

while ~strcmp(ifit_options.line, 'exit') && ~strcmp(ifit_options.line, 'return')
  % handle specific commands (to override limitations from stand-alone)
  if strcmp(ifit_options.line, 'help')        % 'help' command ---------------------------------
    disp('Enter any Matlab/iFit command.');
    disp('      Use ''run script.m'' to execute a script from a file.');
    disp('      Control statements are allowed (for/while loops, tests, switch/case...).');
    disp('Keys: Arrow-Up/Down  Navigate in command history.');
    disp('      Ctrl-C         Exit (same as ''exit'' or ''return'' commands.');
    disp('Help: Type ''doc(iData,''iFit'')'' or see <ifit.mccode.org> <ifit-users@mccode.org>.');
    disp(' ');     
    ifit_options.line = 'doc(iData,''iFit'');';
  elseif strncmp(ifit_options.line,'run ', 4) % 'run' command ----------------------------------
    ifit_options.line = strtrim(ifit_options.line(5:end));
    if ~isempty(dir(ifit_options.line))
      [ifit_options.p,ifit_options.f,ifit_options.e]=fileparts(ifit_options.line);
      if isempty(ifit_options.p), ifit_options.p=pwd; end
      ifit_options.line = fullfile(ifit_options.p, [ifit_options.f ifit_options.e ]);
      disp([ 'Reading file ' ifit_options.line ':' ]);
      ifit_options.fid=fopen(ifit_options.line, 'r');
      ifit_options.line=fread(ifit_options.fid, Inf, 'uint8=>char')';
      fclose(ifit_options.fid);
      ifit_options=rmfield(ifit_options, 'fid');
      ifit_options=rmfield(ifit_options, 'p');
      ifit_options=rmfield(ifit_options, 'f');
      ifit_options=rmfield(ifit_options, 'e');
    else
      disp(['Can not open file ' ifit_options.line]);
    end
  end
  
  % now do the work (evaluate what to do) --------------------------------------
  try
    if ischar(ifit_options.line) && ~isempty(ifit_options.line), eval(ifit_options.line); end
  catch
    if length(ifit_options.line) > 250
      ifit_options.line = [ ifit_options.line(1:250) ' ...' ];
    end
    disp('Error when evaluating:')
    disp(ifit_options.line)
    disp(lasterr)
    ifit_options.line = '';
  end
  
  % collect next command to execute
  if ~isempty(varargin) % from command line ------------------------------------
    ifit_options.not_a_file = 0;
    ifit_options.line = varargin{1}; varargin(1) = [];
    disp([ 'iFit:' num2str(ifit_options.index) '>> argument ' ifit_options.line ]);
    % specific case of imported arguments from the command line
    if (ifit_options.line(1)=='"' && ifit_options.line(end)=='"')
      % a string explicitly indicated as such
      ifit_options.line=ifit_options.line(2:(end-1));
      this{end+1} = ifit_options.line;
    elseif (ifit_options.line(1)=='''' && ifit_options.line(end)=='''')
      % an expression explicitly indicated as such
      ifit_options.line=ifit_options.line(2:(end-1));
      try
        ifit_options.line = eval(ifit_options.line);
      catch
        disp('Error when evaluating argument:')
        disp(ifit_options.line)
        disp(lasterr)
      end
      this{end+1} = ifit_options.line;
    elseif strcmp(ifit_options.line, '--save') || strcmp(ifit_options.line, '-s')
      ifit_options.save=1;
    elseif strncmp(ifit_options.line, '--save=', 7)
      options.save=ifit_options.line(8:end);
    elseif strcmp(ifit_options.line, '--exit') || strcmp(ifit_options.line, '-e')
      ifit_options.exit=1;
    elseif strcmp(ifit_options.line, '--help') || strcmp(ifit_options.line, '-h')
      ifit_options.t = textwrap({ version(iData,'contrib') },80);
      fprintf(1, '%s\n', ifit_options.t{:});
      NL = sprintf('\n');
      ifit_options.t= [ 'Usage:  ifit [options] argruments commands ' NL ...
        NL ...
        '  iFit executable from command line.' NL ...
        '  Project page at <http://ifit.mccode.org>' NL ...
        NL ...
        '  arguments:' NL ...
        '    Any file name, including directories (imports all recursively).' NL ...
        '      Files can also be given as distant URLs (via http and ftp), and as compressed' NL ...
        '      content (zip, tar, gzip).' NL ...
        '      Script files (.m) are executed, other files are imported.' NL ...
        '    Any numerical value, including arrays given as e.g. ''[ 1 2 ...]''.' NL ...
        '    The ''empty'' argument is obtained from the ''[]'' argument' NL ...
        '    Any argument given between double quotes is used as a string argument.' NL ...
        '    Any argument given between simple quotes is evaluated as a an expression.' NL ...
        NL ...
        '  Once all arguments have been processed, they are stored in the ''this''' NL ...
        '    variable (cell array) for further use, and the program enters in interactive' NL ...
        '    mode except when the --exit argument is specified.' NL ...
        '  In interactive mode, any Matlab command can be entered.' NL ...
        NL ...
        '  options:' NL ...
        '  --exit or -e' NL ...
        '      exits immediately after all execution of command line arguments.' NL ...
        '  --save or -s or --save=FILE' NL ...
        '      save the all variables when commands have been executed.' NL ...
        NL ...
        '  Examples:' NL ...
        '    ifit --save file1.*  subplot ' NL NL ];
      disp(ifit_options.t);
      ifit_options=rmfield(ifit_options, 't');
      if isempty(varargin) && ~isempty(this) % last argument has just been processed
        disp('Info: all imported arguments have been stored in ''this''.');
        disp(this);
      end
    elseif ~isempty(str2num(ifit_options.line))
      % numerical value(ifit_options.line) as a matrix
      this{end+1} = str2num(ifit_options.line);
    else
      % check if a file import/conversion to iData is needed
      % apparently not a file name ? -> store string as is
      % must not include filesep, not start with 'file','http','ftp',...
      % and must not exist in the path
      ifit_options.not_a_file =    isempty(strfind(ifit_options.line, filesep)) ...
                   && isempty(strfind(ifit_options.line, '\')) ...
                   && isempty(strfind(ifit_options.line, '/')) ...
                   && isempty(strncmp(ifit_options.line,{'file:','http:'},5)) ...
                   && isempty(strncmp(ifit_options.line,'ftp:',4)) ...
                   && isempty(dir(ifit_options.line));
      if ~ifit_options.not_a_file % probably a file name
        % check if the file is a script
        [ifit_options.p,ifit_options.f,ifit_options.e]=fileparts(ifit_options.line);
        if strcmp(ifit_options.e, 'm')
          % request to execute the script
          ifit_options.line = [ 'run ' ifit_options.line ];
          ifit_options.not_a_file = nan;  % will retain the line as a command
        elseif ismethod(iData, ifit_options.line) || any(exist(ifit_options.line) == [ 3 5 6 ])
          ans = feval(ifit_options.line, this{:})
        else
          this{end+1} = iData(ifit_options.line);
          
        end
        ifit_options=rmfield(ifit_options, 'p');
        ifit_options=rmfield(ifit_options, 'f');
        ifit_options=rmfield(ifit_options, 'e');
      end

    end
    if ~isnan(ifit_options.not_a_file)
      ifit_options.line = '';
    end
    ifit_options=rmfield(ifit_options, 'not_a_file');
    ifit_options.index=ifit_options.index+1;
    
  elseif ifit_options.exit
    ifit_options.line = 'exit';
  else
    if ~isdeployed
      ifit_options.line = input([ 'iFit:' num2str(ifit_options.index) '>> ' ],'s');
    else
      ifit_options.line = input([ 'iFit:' num2str(ifit_options.index) ],'s');
    end
  end
  if ~isempty(ifit_options.line)
    ifit_options.index=ifit_options.index+1;
  end
end

% auto save ?
if isnumeric(ifit_options.save) && ifit_options.save
  save('ifit.mat','-v7.3');
elseif ischar(ifit_options.save)
  save(ifit_options.save, '-v7.3');
end

close all
disp([ '** Ending iFit on ' datestr(now) ])
disp('** Thanks for using iFit <ifit.mccode.org> **');

