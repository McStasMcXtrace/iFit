function filename = ResLibCal_Saveas(filename, EXP, flag)
% ResLibCal_Saveas(filename, EXP): save the GUI configuration to an EXP/ResLib file
%   supports saving in : INI  ResLibCal configuration
%                        DAT  Rescal  Cooper-Nathans legacy (only numbers)
%                        PAR  Rescal5 Cooper-Nathans (with comments)
%                        CFG  Rescal5 Popovici (with comments)
%                        RES  ResCal/ResTrax named list
%                        RTX  ResTrax configuration
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
          '*.dat','ResCal legacy Cooper-Nathans (only numbers, *.dat)' ; ...
          '*.par','ResCal5 Cooper-Nathans (with comments, *.par)' ; ...
          '*.cfg','ResCal5 Popovici configuration (with comments, *.cfg)' ; ...
          '*.res','Cooper-Nathans+Popovici named list (*.res)' ; ...
          '*.rtx','ResTrax legacy (*.rtx)' ; ...
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
      EXP.Rescal.version=5.06; % for ResTrax compatibility
      str = class2str(' ', EXP.ResCal, 'no comment');
      str = [ '% generated automatically on ' datestr(now) ' with ifit.mccode.org ResLibCal' sprintf('\n') ...
        strrep(str, ';', '')] ; % remove trailing ';'
      description = 'named list (Cooper-Nathans+Popovici)';
    elseif any(strcmp(e,'.dat')) && isfield(EXP,'ResCal')
      % Rescal instrument parameters for Cooper-Nathans
      
      % get Rescal labels
      [dummy, labels] = ResLibCal_EXP2RescalPar([]);
      % format the parameter values into a list of numbers in a string, with labels following
      str = struct2cell(EXP.ResCal);
      
      % the legacy Rescal format only saves for Copper-Nathans, i.e. 42 parameters
      if length(str) > 42, str=str(1:42); labels = labels(1:42); end

      str = sprintf('%15.5f\n', str{:});
      description = 'ResCal (Cooper-Nathans) legacy';
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
    elseif any(strcmp(e,'.rtx')) && isfield(EXP,'ResCal')
      ResCal = EXP.ResCal;
      str={ 'Title (max.60 characters): ' ,... 
        [ 'ResLibCal ' datestr(now) ], ...
          'Source (shape,diameter,width,height): ', ...
				[ '0 ' num2str([ sqrt(sum([ResCal.WB ResCal.HB].^2)) ResCal.WB ResCal.HB]) ], ...
					'n-guide (use,distance,length,hor1,hor2,ver1,ver2,ro[m-1],mh,mv[Ni nat.],refh,refv):', ...
					'0 10.0 6300.0 2.5 2.5 15.0 15.0 2.4E-4 1.0 1.0 1.0 1.0', ...
					'Monochromator (chi,aniz.,poiss.,thick.,height,length,segments hor. & vert.):', ...
			  [ '0.0 1.0 .033 ' num2str([ResCal.TM ResCal.HM ResCal.WM]) ' 1 3' ], ...
					'Analyzer (chi,aniz.,poiss.,thick.,height,length,segments hor. & vert.):', ...
			  [ '0.0 1.0 0.3 '  num2str([ResCal.TA ResCal.HA ResCal.WA]) ' 1 3' ], ...
					'Detector (shape,diameter,width,height): ', ...
				[ '1 0. ' num2str([ResCal.WD ResCal.HD]) ' 3.0 5.0' ], ...
					'Distances (L1,L2,L3,L4):', ...
					num2str([ ResCal.L1 ResCal.L2 ResCal.L3 ResCal.L4 ]), ...
					'1st collimator (distance,length,hor1,hor2,ver1,ver2,ro[m-1],mh,mv[Ni nat.],refh,refv):', ...
				[ num2str([ResCal.L1/4 ResCal.L1/2]) ' 8.05 5.0 9.05 11.0 0.0 0.0 0.0 1.0 1.0' ], ...
					'2nd collimator (distance,length,hor1,hor2,ver1,ver2,ro[m-1],mh,mv[Ni nat.],refh,refv):', ...
				[ num2str([ResCal.L1/4 ResCal.L1/2]) ' 4.0 4.0 7.0 7.0 0.0 0.0 0.0 1.0 1.0' ], ...
					'3rd collimator (distance,length,hor1,hor2,ver1,ver2,ro[m-1],mh,mv[Ni nat.],refh,refv):', ...
				[ num2str([ResCal.L1/4 ResCal.L1/2]) ' 4.0 4.0 7.0 7.0 0.0 0.0 0.0 1.0 1.0' ], ...
					'4th collimator (distance,length,hor1,hor2,ver1,ver2,ro[m-1],mh,mv[Ni nat.],refh,refv):', ...
				[ num2str([ResCal.L1/4 ResCal.L1/2]) ' 4.0 4.0 12.0 12.0 0.0 0.0 0.0 1.0 1.0' ], ...
					'COLLIMATORS  1  1  1  1' };
		  str = sprintf('%s\n', str{:});
			description = 'ResTrax legacy';
			disp('WARNING: the generated ResTrax file does not contain the full configuration');
    else
      NL = sprintf('\n');
      str = [ '% ResLibCal configuration script file ' NL ...
            '%' NL ...
            '% Matlab ' version ' m-file ' filename NL ...
            '% generated automatically on ' datestr(now) ' with ifit.mccode.org ResLibCal' NL...
            '% The Resolution function is indicated as the "resolution.RMS" field ' NL ...
            '% (in lattice rlu), "resolution.RM" is in [Angs-1] in [QxQyQzE].' NL ...
            '% If you ever edit this file manually, please modify "config" rather than "config.ResCal".' NL ...
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
