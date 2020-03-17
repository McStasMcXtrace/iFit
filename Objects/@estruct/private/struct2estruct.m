function b=struct2estruct(a, varargin)
% struct2estruct converts a structure into an estruct
%   STRUCT2ESTRUCT(struct) a new estruct object is created, which contains
%   the initial struct
%
%   STRUCT2ESTRUCT(struct, org) when a 2nd argument is given as an estruct, it 
%   is used as original object, and updated with the structure
  persistent fb

  if isempty(fb), fb=fieldnames(estruct); end
  
  if ~isstruct(a), b=[]; return; end
  if numel(a) > 1
    b = [];
    for index=1:numel(a)
      b = [ b ; struct2estruct(a(index), varargin{:}) ];
    end
    return
  end

  f  = fieldnames(a);
  if nargin == 1
    b = estruct;
  else
    b = varargin{1};
  end
  if isfield(a, 'Data')   % start by storing the raw Data
    b.Data = a.Data;
  end
  for index=1:length(f) % store the estruct static fields
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
  
  if isfield(a, 'Command')
    b.Command = a.Command;
  end
  % transfer some standard fields from iLoad when empty/non-existent
  for f={'Source','Date','Label','Format','User','Loader'}
    if isfield(a, f{1}) && ~isempty(a.(f{1})) ...
    && (~isfield(b, f{1}) || (isempty(getalias(b,f{1})) && isempty(getalias(b,f{1})))) ...
      if isfield(b, f{1})
        b.(f{1}) = a.(f{1});
      else
        set(b, f{1}, a.(f{1}));
      end
    end
  end
  if isfield(a, 'Title') && ~isempty(a.Title), b.Name = a.Title; end
  if isfield(a, 'Loader') && isfield(a.Loader, 'postprocess') && ~isempty(a.Loader.postprocess)
    set(b, 'postprocess', a.Loader.postprocess);
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

