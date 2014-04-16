function filename = ResLibCal_Saveas(filename, EXP, flag)
% ResLibCal_Saveas(filename, EXP): save the GUI configuration to an EXP/ResLib file
%   supports saving in : INI  ResLibCal configuration
%                        PAR  Rescal Cooper-Nathans
%                        CFG  Rescal5 Popovici (also needs PAR)
%                        RES  named list
%
%   CFG file will also save a PAR file, and vice versa, except when a 3rd argument is specified.
% Return:
%  filename: name of file written

% Calls: ResLibCal_Compute, class2str

  if nargin < 1, filename =''; end
  if nargin < 2, EXP=''; end

  if isempty(filename)
    [filename, pathname] = uiputfile( ...
         {'*.m',  'M-files (*.m)'; ...
          '*.ini','ResLibCal configuration (*.ini)' ; ...
          '*.par','ResCal Cooper-Nathans configuration (*.par)' ; ...
          '*.cfg','ResCal Popovici configuration (*.cfg)' ; ...
          '*.res','ResCal named list (*.res)' ; ...
          '*.*',  'All Files (*.*)'}, ...
          'Save ResLibCal configuration as ...');
    if isempty(filename) || all(filename == 0), return; end
    filename = fullfile(pathname, filename);
  end

  if isempty(EXP)
    EXP = ResLibCal_Compute; % collect current settings and compute
  end
  if ~isempty(EXP)
    % do we reduce output to ResCal only ? based on extension match
    [p,f,e] = fileparts(filename);
    if any(strcmp(e,'.res')) && isfield(EXP,'ResCal')
      % 'NAME = value' list file
      str = class2str(' ', EXP.ResCal, 'no comment');
      str = [ '% generated automatically on ' datestr(now) ' with ifit.mccode.org ResLibCal' sprintf('\n') ...
        strrep(str, ';', '')] ; % remove trailing ';'
      description = 'named list (Cooper-Nathans+Popovici)';
    elseif any(strcmp(e,'.par')) && isfield(EXP,'ResCal')
      % Rescal instrument parameters for Cooper-Nathans
      
      % get Rescal labels
      [dummy, labels] = ResLibCal_EXP2RescalPar([]);
      % format the parameter values into a list of numbers in a string, with labels following
      str = struct2cell(EXP.ResCal);
      
      % the legacy Rescal format only saves for Copper-Nathans, i.e. 42 parameters
      if length(str) > 42, str=str(1:42); labels = labels(1:42); end
      
      % assemble values and labels in a unique cell
      s=cell(1,2*length(str)); s(1:2:end) = str; s(2:2:end) = labels;
      str = [ sprintf('%-15g %% %s\n', s{:}) ...
         '% ResCal parameters for Cooper-Nathans. ' datestr(now) ' ifit.mccode.org ResLibCal' sprintf('\n') ];
      description = 'ResCal (Cooper-Nathans)';
    elseif any(strcmp(e,'.cfg')) && isfield(EXP,'ResCal')
      % Rescal instrument parameters for Popovici
      
      % get Rescal labels
      [dummy, labels] = ResLibCal_EXP2RescalPar([]);
      % format the parameter values into a list of numbers in a string, with labels following
      str = struct2cell(EXP.ResCal);
      
      % the Rescal5 'cfg' file contains the Popovici additional 27 parameters
      if length(str) > 42, str=str(43:69); labels = labels(43:69); end
      
      % assemble values and labels in a unique cell
      s=cell(1,2*length(str)); s(1:2:end) = str; s(2:2:end) = labels;
      str = [ sprintf('%-15g %% %s\n', s{:}) ...
         '% ResCal parameters for Popovici. ' datestr(now) ' ifit.mccode.org ResLibCal' sprintf('\n') ];
      
      description = 'ResCal5 (Popovici)';
    else
      str = [ '% ResLibCal configuration script file ' sprintf('\n') ...
            '%' sprintf('\n') ...
            '% Matlab ' version ' m-file ' filename sprintf('\n') ...
            '% generated automatically on ' datestr(now) ' with ifit.mccode.org ResLibCal' sprintf('\n') ...
            '% The Resolution function is indicated as the "resolution.RMS" field ' sprintf('\n') ...
            '% (in lattice rlu), "resolution.RM" is in [Angs-1] in [QxQyQzE].' sprintf('\n') ...
            '% If you ever edit this file manually, please modify "config" rather than "config.ResCal".' sprintf('\n') ...
            class2str('config', EXP) ];
      description = 'ResLibCal (Cooper-Nathans+Popovici)';
    end
    
    [fid, message]=fopen(filename,'w+');
    if fid == -1
      warning(['Error opening file ' filename ' to save ' description ' configuration.' ]);
      filename = [];
    else
      fprintf(fid, '%s', str);
      fclose(fid);
      disp([ '% Saved ' description ' configuration into ' filename ]);
      
      % for .par and .cfg files, we must also save the other complementary file
      if nargin < 3
        if any(strcmp(e,'.cfg'))
          % CFG -> save PAR
          ResLibCal_Saveas(fullfile(p, [f '.par' ]), EXP, 'no recursive!');
        elseif any(strcmp(e,'.par'))
          % PAR -> save CFG
          ResLibCal_Saveas(fullfile(p, [f '.cfg' ]), EXP, 'no recursive!');
        end
      end
    end
  end
% end ResLibCal_Saveas
