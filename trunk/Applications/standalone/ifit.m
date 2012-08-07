function ifit(varargin)
% Usage:  ifit [options] argruments commands 
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
%  To execute multiple lines and control statements (if/while/for...), write them
%    in a script and type 'run <script>'.
%
%  options:
%  --exit or -e
%      exits immediately after all execution of command line arguments.
%  --save or -s or --save=FILE
%      save the all variables when commands have been executed.
%  --run=SCRIPT or -r=SCRIPT
%      executes the SCRIPT when starting.
%
%  Examples:
%    ifit --save file1.*  subplot 

% Manual build:
% Better use the 'make.sh' script, or
% addpath(genpath('/home/farhi/svn/Matlab/iFit/trunk'))
% mcc -m ifit -a /home/farhi/svn/Matlab/iFit/trunk
% buildmcr('.')

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
disp('Type ''help'' to learn how to use this software. Type ''exit'' or Ctrl-C to exit.');
if ispc
  disp('WARNING: under Windows platforms, file names containing spaces, such as "My Documents" ')
  disp('         are not well supported. Rename files and move directories to other locations.')
end
disp(' ')

ifit_options.line     ='';     % the current line to execute
ifit_options.index    =1;      % the index of the input
this                  ={};     % the buffer from the command line
ifit_options.save     =0;
ifit_options.exit     =0;

while ~strcmp(ifit_options.line, 'exit') && ~strcmp(ifit_options.line, 'return')
  ifit_options.line = strtrim(ifit_options.line);
  % handle specific commands (to override limitations from stand-alone)
  if strcmp(strtok(ifit_options.line), 'doc')
    ifit_options.line = [ 'help' ifit_options.line(4:end) ];
  end
  if strcmp(strtok(ifit_options.line), 'help')        % 'help' command ---------------------------------
    if length(ifit_options.line) > 4
    disp(ifit_options.line)
      [ifit_options.t, ifit_options.line] = strtok(ifit_options.line);
      ifit_options=rmfield(ifit_options, 't');
      ifit_options.line      = strtrim(ifit_options.line);
      if ~isempty(ifit_options.line), 
        web(ifit_options.line); 
        ifit_options.line = '';
      end
    else
      disp('Enter any Matlab/iFit command.');
      disp('      Use ''run script.m'' to execute a script from a file.');
      disp('      Control statements are allowed (for/while loops, tests, ');
      disp('      switch/case...) when they span on one line, or in scripts.');
      disp('Keys: Arrow-Up/Down  Navigate in command history.');
      disp('      Ctrl-C         Exit (same as ''exit'' or ''return'' commands.');
      disp('Help: Type ''doc(iData,''iFit'')'' ');
      disp('      see or <ifit.mccode.org> <ifit-users@mccode.org>.');
      disp('To import some data, use e.g. d=iData(''filename'');');
      disp('To create a model, use e.g. f=iFunc(''expression''); ');
      disp('  or type fits(iFunc) to get a list of available models.');
      disp('To fit a model to data, use e.g. fits(f,d)');
      disp('Data and Models can be manipulated (+-/*...) using the Matlab syntax.');
      disp(' ');     
      disp('Matlab is a registered trademark of The Mathworks Inc.');
      disp('Source code for this software is available at <ifit.mccode.org>.')
      disp('Matlab help is fully available at <http://www.mathworks.com/help/techdoc>.');
      ifit_options.line = 'doc(iData,''iFit''); disp('' '');';
    end
  elseif strncmp(ifit_options.line,'run ', 4) % 'run' command ----------------------------------
    ifit_options.line = strtrim(ifit_options.line(5:end));
    if ~exist(ifit_options.line) && exist([ ifit_options.line '.m' ])
      ifit_options.line = [ ifit_options.line '.m' ];
    end
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
      disp(['Error: iFit: Can not open file ' ifit_options.line]);
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
      ifit_options.save='ifit.mat';
    elseif strncmp(ifit_options.line, '--save=', 7)
      options.save=ifit_options.line(8:end);
    elseif strncmp(ifit_options.line, '--run=', 6)
      options.save=[ 'run ' ifit_options.line(7:end) ];
    elseif strncmp(ifit_options.line, '-r=', 3)
      options.save=[ 'run ' ifit_options.line(4:end) ];
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
        '  --run=SCRIPT or -r=SCRIPT' NL ...
        '      executes the SCRIPT when starting.' NL ...
        NL ...
        '  Examples:' NL ...
        '    ifit --save file1.*  subplot ' NL NL ];
      disp(ifit_options.t);
      ifit_options=rmfield(ifit_options, 't');
      if isempty(varargin) && ~isempty(this) % last argument has just been processed
        disp('Info: all imported arguments have been stored in cell array ''this''.');
        disp('      access them with e.g. this{1} ...');
        disp(this);
        clear varargin
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
          ans = this{end}
          if isempty(varargin) && length(this) <= 20 % last file imported. Plot all imported data sets
            subplot(this{cellfun('isclass',this,'iData')});
          end
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
    ifit_options.line = input([ 'iFit:' num2str(ifit_options.index) ' ' ],'s');
  end
  if ~isempty(ifit_options.line)
    ifit_options.index=ifit_options.index+1;
  end
end

% auto save ?
if ischar(ifit_options.save)
  save(ifit_options.save, '-v7.3');
end

close all
disp([ '** Ending iFit on ' datestr(now) ])
disp('** Thanks for using iFit <ifit.mccode.org> **');

