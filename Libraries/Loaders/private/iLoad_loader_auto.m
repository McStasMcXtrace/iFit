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

  loaders={};
  
  [loaders_ext, formats] = iLoad_loader_extension(formats, ext);
  [loaders_pat, formats] = iLoad_loader_patterns(formats, file_start);
  [loaders_txt, formats] = iLoad_loader_text(formats);

  if isbinary
  % binary header: loaders = [ extension, pattern matches, pattern char empty ]
    loaders = [ loaders_ext loaders_pat formats ];
  else
  % text header:   loaders = [ pattern matches, extension, pattern char empty ]
    loaders = [ loaders_pat loaders_ext loaders_txt formats ];
  end
  
end % iLoad_loader_auto

% ------------------------------------------------------------------------------
function [loaders, others] = iLoad_loader_extension(formats, ext)
% iLoad_loader_extension: identify loaders for which any registered extension match 'ext'
  exts = cellfun(@(c)getfield(c, 'extension'), formats, 'UniformOutput',false);
  loaders={}; others={};
  for index=1:numel(exts)
    if any(strcmp(ext, exts{index})) 
        loaders{end+1} = formats{index}; 
    else others{end+1} = formats{index}; end
  end
  
end % iLoad_loader_extension

function [loaders, others] = iLoad_loader_patterns(formats, file_start)
% iLoad_loader_patterns: identify loaders for which all patterns appear in file_start 
  pats = cellfun(@(c)getfield(c, 'patterns'), formats, 'UniformOutput',false);
  loaders={}; others={};
  for index=1:numel(pats)
    if ~iscellstr(pats{index}) && ~ischar(pats{index}), continue; end
    found=true;
    for pat=cellstr(pats{index})
      match = strfind(file_start, pat{1});
      if isempty(match), found=false; break; end
    end
    if found
        loaders{end+1} = formats{index};
    else others{end+1} = formats{index}; end
  end
  
end % iLoad_loader_patterns

function [loaders, others] = iLoad_loader_text(formats)
% iLoad_loader_text: identify loaders for which have empty char patterns (text readers)
  istext  = cellfun(@(c)getfield(c, 'istext'), formats);
  pats    = cellfun(@(c)getfield(c, 'patterns'), formats, 'UniformOutput',false);
  emptypat= cellfun(@(c)and(~isempty(c),isempty(char(c))), pats);
  istext = istext & emptypat;
  loaders = formats(find(istext));
  others  = formats(find(~istext));
end % iLoad_loader_text

