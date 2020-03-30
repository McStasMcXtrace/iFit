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

  % restrict loaders using format extensions, also retain formats without extension
  %   this filter is needed to restrict formats in a fast way. 
  [loaders,selected] = iLoad_loader_extension(formats, ext);

  % restrict loaders that match patterns (or no patterns). 
  % Keep extension exact match on top
  loaders = iLoad_loader_patterns(formats, file_start, isbinary, ext, selected);

end % iLoad_loader_auto

% ------------------------------------------------------------------------------
function [loaders,selected] = iLoad_loader_extension(formats, ext)
% iLoad_loader_extension: identify loaders for which any registered extension match 'ext'
%   or loader has no extension (all match), or '*'

  persistent exts
  if isempty(exts)
    exts = cellfun(@(c)getfield(c, 'extension'), formats, 'UniformOutput',false);
  end
  
  done = zeros(size(exts));
  loaders={};
  % first we add formats that exactly match extensions
  for index=1:numel(exts)
    if ~isempty(exts{index}) && isempty(ext), found = false; 
    else found=any(strcmpi(ext, exts{index}));
      if found, loaders{end+1} = formats{index}; done(index)=1; end
    end
  end
  % then we add other ones
  for index=1:numel(exts)
    if ~done(index)
      if ~isempty(exts{index}) && isempty(ext), found = false;
      else found=isempty(exts{index}) || (iscellstr(exts{index}) && any(strcmp(exts{index},'*')));
        if found, loaders{end+1} = formats{index}; done(index) = 1; end
      end
      
    end
  end
  selected = find(done);
  
end % iLoad_loader_extension

% ------------------------------------------------------------------------------
function loaders = iLoad_loader_patterns(formats, file_start, isbinary, ext, selected)
% iLoad_loader_patterns: identify loaders for which all patterns appear in file_start 

  % keep when patterns found in file_start
  % or format has patterns=='' and is text
  % or no pattern at all
  persistent exts0 pats0 istext0
  if isempty(exts0)
    f = iLoad_config_load;
    exts0   = cellfun(@(c)getfield(c, 'extension'), f.loaders, 'UniformOutput',false);
    pats0   = cellfun(@(c)getfield(c, 'patterns'),  f.loaders, 'UniformOutput',false);
    istext0 = cellfun(@(c)getfield(c, 'istext'),    f.loaders);
  end
  % restrict formats to selected set (from extension)
  exts  = exts0(selected);
  pats  = pats0(selected);
  istext= istext0(selected);
  formats=formats(selected);
  
  done   = zeros(size(istext));
  nopat  = zeros(size(istext));
  loaders= {};
  % first put formats that match patterns
  for index=1:numel(pats)
    found = false;
    if (iscellstr(pats{index}) || ischar(pats{index})) && ~isbinary
      found=true;
      for pat=cellstr(pats{index})
        match = regexp(file_start, pat{1});
        if isempty(match), found=false; nopat(index)=1; break; end % a pattern was not found
      end
      if found, 
        loaders{end+1} = formats{index};
        done(index)=1;
      end
    end
  end
  % we keep formats that exactly match extension (but not if patterns do not match)
  for index=1:numel(exts)
    if ~done(index) && ~nopat(index)
      if ~isempty(exts{index}) && isempty(ext), found = false; 
      else found=any(strcmpi(ext, exts{index})); end
      if found
        loaders{end+1} = formats{index}; done(index)=1; 
      end
    end
  end
  % then put loaders without specific pattern
  for index=1:numel(pats)
    if ~isbinary && ~done(index)
      found = false;
      if (iscellstr(pats{index}) || ischar(pats{index})) && all(strcmp(pats{index},''))
        found = istext(index);  % text format without specific pattern
        if found, 
          loaders{end+1} = formats{index};
          done(index)=1;
        end
      end
      
    end
  end
  % then put loaders without patterns at all (including binary ones)
  for index=1:numel(pats)
    if isempty(pats{index}) && ~done(index), loaders{end+1} = formats{index}; end
  end

end % iLoad_loader_patterns

