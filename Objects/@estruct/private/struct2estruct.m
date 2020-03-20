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

  if nargin == 1
    b = estruct;
  elseif isa(varargin{1}, 'estruct')
    b = varargin{1};
  end
  
  if isfield(a, 'Title') && ~isfield(a, 'Name')
    a.Name = a.Title; a = rmfield(a, 'Title');
  end

  % transfer the fields (except protected ones)
  for f=fieldnames(a)'
    if any(strcmp(f{1}, b.properties_Protected)), continue; end % ignore protected
    if ismethod(b, f{1})
      warning([ mfilename ': skipping input struct field ' f{1} ' as it matches a method name' ])
    elseif isfield(b, f{1})
      b.(f{1}) = a.(f{1});
    else
      set(b, f{1}, a.(f{1})); % add new property
    end
  end
  
  % handle some more specific stuff
  % Attributes Headers MetaData Loader postprocess
  if isfield(a, 'Headers')
    label(b, 'Headers', 'Headers (text)' );
  end
  if isfield(a, 'Attributes')
    label(b, 'Attributes', 'Attributes (text)' );
  end
  if isfield(b.Data,'MetaData') && ~isfield(b, 'MetaData')
    set(b, 'MetaData', 'Data.MetaData');
  end
  if isfield(b, 'MetaData')
    [pathname,filename,ext] = fileparts(b.Source);
    label(b, 'MetaData', [ 'MetaData from ' filename ext ] );
    try
      b=load_clean_metadata(b);
    end
  end
  if ~isfield(b,        'postprocess') &&  isfield(a, 'Loader') ...
  &&  isfield(a.Loader, 'postprocess') && ~isempty(a.Loader.postprocess)
    set(b, 'postprocess', a.Loader.postprocess);
  end
  
  if ~iscellstr(b.Command)
    b.Command = { b.Command };
  end
  
  if isempty(b.Command), b.Command= cellstr('estruct(<struct>)'); end
  
  % ------------------------------------------------------------------------------
  
function a=load_clean_metadata(a)
% LOAD_CLEAN_METADATA test each field of MetaData and search for equal aliases
  this = a.MetaData;
  meta_names = fieldnames(this);
  alias_names= getalias(a);
  % treat each MetaData
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
  a.MetaData = this;

