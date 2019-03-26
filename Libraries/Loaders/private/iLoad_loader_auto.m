function [loaders, isbinary] = iLoad_loader_auto(file, config)
% function to determine which parser to use to analyze content

  verbose = 0;  % set it to 1 to display Loader auto-selection from extension/pattern

  % get the configuration
  % config  = iLoad('','load config');
  
  loaders = config.loaders;
  isbinary= 0;
  % read start of file
  [fid,message] = fopen(file, 'r');
  if fid == -1
    fprintf(1, 'iLoad: %s: %s. Check existence/permissions.\n', file, message );
    error([ 'Could not open file ' file ' for reading. ' message '. Check existence/permissions.' ]);
  end
  file_start = fread(fid, 1000, 'uint8=>char')';
  
  % loop to test each format for patterns
  formats = loaders;
  loaders ={}; % will hold those formats which are applicable
  loaders_count=0;
  % identify by extensions
  [dummy, dummy, fext] = fileparts(file);
  if ~isempty(fext) && fext(1)=='.', fext(1)=[]; end
  % check if this is a text file
  if length(find(file_start >= 32 & file_start < 127))/length(file_start) < 0.4
    isbinary = 1; % less than 90% of printable characters
  else
    % clean file start with spaces, remove EOL
    file_start = [ file_start fread(fid, 9000, 'uint8=>char')' ];
    file_start(isspace(file_start)) = ' ';
  end
  fclose(fid);
  
  % Selection procedure:
  %   when pattern is empty char (text) restrict to 'istext' loaders
  
  %   when patterns exist, not empty: find them -> select loader
  
  %   when extension matches: -> select loader
  
  
  % identify by patterns
  for index=1:length(formats)
    loader = formats{index};
    if ~isstruct(loader), continue; end

    if ~isfield(loader,'patterns'), loader.patterns=[]; end
    
    patterns_found = 0;
    if verbose, fprintf(1,'iLoad: method %s: file %s: analyzing\n', loader.name, file); end
    
    % the loader is selected if: 
    %   loader has extension and extension matches, no patterns to search
    if ~patterns_found && (...
      (isempty(loader.patterns) && ~ischar(loader.patterns) && ~iscellstr(loader.patterns)) ...
      || isnumeric(loader.patterns))
      if ~isfield(loader,'extension'), ext=''; 
      else ext=loader.extension; end
      if ischar(ext) && ~isempty(ext), ext={ ext }; end
      if ~isempty(ext) && ~isempty(fext) 
        if any(strcmpi(fext, ext))
          patterns_found  = 1;  % extension does match
          if verbose, disp([ 'iLoad: method ' loader.name ': ' file ': extension ' ext{1} ' matches' ]); end
        end
      end
    end
      
    %   loader has patterns and not binary and patterns match file_start, 
    %   whatever be the extension
    % ~isbinary may be removed in case it suppresses e.g. text/binary formats such as EDF, ADSC, ...
    if ~isempty(loader.patterns) && ischar(loader.patterns)
      loader.patterns=cellstr(loader.patterns); 
    end
    
    if ~patterns_found && ~isempty(loader.patterns) && iscellstr(loader.patterns) && ~isbinary 
    
      % check all patterns in text file
      all_match = 1;
      for index_pat=1:numel(loader.patterns) % all patterns must match
        if isempty(regexpi(file_start, loader.patterns{index_pat}, 'once'))
          all_match=0;     % at least one pattern does not match
          if verbose, fprintf(1,'iLoad: method %s: file %s: at least one pattern does not match (%s)\n', loader.name, file, loader.patterns{index_pat}); end
          break;
        end
      end % for patterns
      if all_match, patterns_found=1; 
        if verbose, disp([ 'iLoad: method ' loader.name ': ' file ': patterns match' ]); end
      end
    end
      
    %   loader has no extension and no patterns
    if ~patterns_found && isempty(fext) && isempty(loader.patterns)
      patterns_found = 1; % we will try all non selective loaders
      if verbose, disp([ 'iLoad: method ' loader.name ': ' file ': no patterns nor extension: select by default' ]); end
    end
      
    if patterns_found
      loaders_count = loaders_count+1;
      loaders{loaders_count} = loader;
    end

  end % for index

end % iLoad_loader_auto
