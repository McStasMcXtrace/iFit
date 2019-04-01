function [loaders, isbinary] = iLoad_loader_auto(file, config)
% function to determine which parser to use to analyze content
  
  formats = config.loaders;
  
  % get extension 
  [~, ~, ext] = fileparts(file);
  if ~isempty(ext) && ext(1)=='.', ext(1)=[]; end
  
  % read start of file
  [fid,message] = fopen(file, 'r');
  if fid == -1
    fprintf(1, 'iLoad: ERROR: %s: %s. Check existence/permissions.\n', file, message );
    error([ 'Could not open file ' file ' for reading. ' message '. Check existence/permissions.' ]);
  end

  % get header (file start) and 'isbinary'
  file_start = fread(fid, 1000, 'uint8=>char')';
  if length(find(file_start >= 32 & file_start < 127))/length(file_start) < 0.4
    isbinary = 1; % less than 40% of printable characters
  else
    % clean file start with spaces, remove EOL
    file_start = [ file_start fread(fid, 9000, 'uint8=>char')' ];
    file_start(isspace(file_start)) = ' ';
    isbinary= 0;
  end
  fclose(fid);

  loaders=formats; % already sorted by decreasing specificity
 
 
  % restrict loaders using format extensions, also retain formats without extension
  loaders = iLoad_loader_extension(loaders, ext);

  % restrict loaders that match patterns (or no patterns)
  loaders = iLoad_loader_patterns(loaders, file_start, isbinary);

end % iLoad_loader_auto

% ------------------------------------------------------------------------------
function loaders = iLoad_loader_extension(formats, ext)
% iLoad_loader_extension: identify loaders for which any registered extension match 'ext'
%   or loader has no extension (all match), or '*'
  exts = cellfun(@(c)getfield(c, 'extension'), formats, 'UniformOutput',false);
  done = zeros(size(exts));
  loaders={};
  % first we add formats that exactly match extensions
  for index=1:numel(exts)
    found=any(strcmp(ext, exts{index}));
    if ~isempty(exts{index}) && isempty(ext), found = false; end
    if found, loaders{end+1} = formats{index}; done(index)=1; end
  end
  % then we add other ones
  for index=1:numel(exts)
    found=isempty(exts{index}) || (ischar(exts{index}) && strcmp(exts{index},'*'));
    if ~isempty(exts{index}) && isempty(ext), found = false; end
    if found && ~done(index), loaders{end+1} = formats{index}; end
  end
  
end % iLoad_loader_extension

function [loaders, others] = iLoad_loader_patterns(formats, file_start, isbinary)
% iLoad_loader_patterns: identify loaders for which all patterns appear in file_start 

  % keep when patterns found in file_start
  % or format has patterns=='' and is text
  % or no pattern at all
  pats   = cellfun(@(c)getfield(c, 'patterns'), formats, 'UniformOutput',false);
  istext = cellfun(@(c)getfield(c, 'istext'),   formats);
  done   = zeros(size(istext));
  loaders={};
  % first put formats that match patterns
  for index=1:numel(pats)
    found = false;
    if (iscellstr(pats{index}) || ischar(pats{index})) && ~isbinary
      found=true;
      for pat=cellstr(pats{index})
        match = regexp(file_start, pat{1});
        if isempty(match), found=false; break; end % a pattern was not found
      end
    end
    if found, 
      loaders{end+1} = formats{index};
      done(index)=1;
    end
  end
  % then put loaders without specific pattern
  for index=1:numel(pats)
    found = false;
    if (iscellstr(pats{index}) || ischar(pats{index})) && ~isbinary && ~done(index)
      if all(strcmp(pats{index},''))
        found = istext(index);  % text format without specific pattern
      end
    end
    if found, 
      loaders{end+1} = formats{index};
      done(index)=1;
    end
  end
  % then put loaders without paterns at all (including binary ones)
  for index=1:numel(pats)
    if isempty(pats{index}) && ~done(index), loaders{end+1} = formats{index}; end
  end
  
end % iLoad_loader_patterns

