function v = get(s, varargin)
% GET    Get object properties.
%  V = GET(S,'PropertyName') 
%    Get the value of the specified property/field in the structure.
%    The object S can be an array.
%
%    The 'PropertyName' can be a full structure path, such as 'field1.field2' in
%    in which case the returned value travelled through.
%
%    When the target property is itself a valid structure path (char), it is also 
%    travelled through before getting the value (see below).
%
%  V = GET(S,'PropertyName1','PropertyName2',...)
%    Get multiple properties.
%
%  V = GET(S, ...,'alias')
%    When the PropertyName points to a string value, it is fetched without 
%    travelling through it (get an alias/link). 
%    This syntax is equivalent to GETALIAS(S, 'PropertyName1',...)
%    In this case, the retrieval allows to access internal or external links, 
%    as well as evaluated expression, with the syntax cases for the 'Value':
%     'field'                           a simple link to an other property 'field'
%     'field1.field2...'                a nested link to an other property
%     'file://some_file_path'           a local file URL
%     'http://some_distant_resource'    an HTTP URL (proxy settings may need to be set)
%     'https://some_distant_resource'   an HTTPS URL (proxy settings may need to be set)
%     'ftp://some_distant_resource'     an FTP URL (proxy settings may need to be set)
%     'matlab: some_expression'         some code to evaluate. 'this' refers to the 
%                                       object itself e.g. 'matlab: this.Signal*2'
%
%    File and URL can refer to compressed resources (zip, gz, tar, Z) which are 
%    extracted on-the-fly. In case the URL/file resource contains 'sections' a search token
%    can be specified with syntax such as 'file://filename#token'.
% 
%  GET(S) displays all object properties.
%
% Example: s=estruct; set(s, 'a', 1); set(s, 'b.c','a','alias'); get(s, 'b.c') == 1 && strcmp(get(s, 'b.c','alias'),'a')
% Version: $Date$ $Version$ $Author$
% See also estruct, fieldnames, findfield, isfield, set, get, getalias, setalias, 
%   getaxis, setaxis

% NOTE: the rationale here is to implement all the logic in subsref and just call it.

  if ~isa(s, 'estruct')
    v = builtin('get', s, varargin{:});
    return
  end

  if nargin == 1, v = fieldnames(s); return; end
  
  % handle array of struct
  v = {};
  if numel(s) > 1
    for index=1:numel(s)
      v{end+1} = get(s(index), varargin{:});
    end
    v = reshape(v, size(s));
    return
  end
  
  % check last argument as 'alias' ? will not follow links, but get aliases directly.
  if ischar(varargin{end}) && any(strcmp(varargin{end}, {'link','alias'}))
       follow = false;  % we get aliases and do not link for get/subsref.
       varargin(end) = [];
  else follow = true;   % we travel through links
  end
  
  % handle array/cell of fields
  for index=1:numel(varargin) % loop on requested properties
    name = varargin{index};
    if ~ischar(name) && ~iscellstr(name)
      error([ mfilename ': GET works with property name argument (char/cellstr). The ' num2str(index) '-th argument is of type ' class(name) ]);
    end
    name = cellstr(name);
    for n_index=1:numel(name)
      if follow
        v{end+1} = subsref(s, struct('type','.', 'subs',name{n_index}));
      else
        v{end+1} = subsref(s, struct('type','()','subs',name{n_index}));
        if ~ischar(v{end}), v{end} = []; end % result must be an alias or []
      end
    end
  end
  if numel(v) == 1, v = v{1}; end
  
