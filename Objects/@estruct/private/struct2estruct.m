function b=struct2estruct(a)
% struct2estruct: converts a structure into an estruct

  persistent fb

  if isempty(fb), fb=fieldnames(estruct); end
  
  if ~isstruct(a), b=[]; return; end
  if numel(a) > 1
    b = [];
    for index=1:numel(a)
      b = [ b ; estruct_struct2estruct(a(index)) ];
    end
    return
  end

  f  = fieldnames(a);
  b  = estruct; 
  if isfield(a, 'Data')   % start by storing the raw Data
    b.Data = a.Data;
  end
  for index=1:length(f)
    if any(strcmp(f{index},fb)) && ~strcmp(f{index}, 'Data')
      b = set(b,f{index}, a.(f{index}));
    end
  end
    
  if ~isfield(a, 'Data')   % store whole file content if possible.
    b.Data = a;
%  else
%    disp(['estruct: warning: could not import all fields from structure.' ]);
  elseif isfield(a, 'Headers')
    b.Data.Attributes = a.Headers;
    b=setalias(b, 'Attributes', 'Data.Attributes', 'Headers (text)' );
   elseif isfield(a, 'Attributes')
    b.Data.Attributes = a.Attributes;
    b=setalias(b, 'Attributes', 'Data.Attributes', 'Headers (text)' );
  end
  if isfield(a, 'Format')
    setalias(b, 'Format', a.Format);
  end
  if isfield(a, 'Command')
    b.Command = a.Command;
  end
  if ~iscellstr(b.Command)
    b.Command = { b.Command };
  end
  
  if isempty(b.Command), b.Command= cellstr('estruct(<struct>)'); end
  try
      [pathname,filename,ext] = fileparts(b.Source);
      if isfield(b.Data, 'MetaData')
        b=setalias(b, 'MetaData', 'Data.MetaData', [ 'MetaData from ' filename ext ]);
        b=load_clean_metadata(b);
      end
  end
  
  % ------------------------------------------------------------------------------
  
function a=load_clean_metadata(a, loaders, filenames)
% test each field of MetaData and search for equal aliases
  this = a.Data.MetaData;
  meta_names = fieldnames(this);
  alias_names=getalias(a);
  %treat each MetaData
  for index=1:length(meta_names)
    if numel(getfield(this, meta_names{index})) > 1000
    for index_alias=1:length(alias_names)
      % is it big and equal to an alias value ?
      if isequal(getfield(this, meta_names{index}), get(a, alias_names{index_alias}))
        % yes: store alias in place of MetaData
        this = setfield(this, meta_names{index}, getalias(a, alias_names{index_alias}));
        break
      end
    end % for
    end % if
  end
  a.Data.MetaData = this;

