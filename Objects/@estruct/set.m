function s = set(s, varargin)
% SET    Set estruct properties.
%
%  V = SET(S,'PropertyName','Value') 
%    Set the value of the specified property/field in the structure.
%    The object S can be an array.
%
%    The 'PropertyName' can be a full structure path, such as 'field1.field2' in
%    in which case the value assigment is made recursive (travel through).
%
%    When the target property is itself a valid structure path (char), it is also 
%    travelled through before assigment (see below).
%
%  V = SET(S,'PropertyName1','Value1','PropertyName2','Value2',...)
%    Set multiple properties.
%
%  V = SET(S, ...,'alias')
%    When the PropertyName points to a string value, it is assigned without 
%    travelling through it (set as alias/link). 
%    This syntax is equivalent to SETALIAS(S, 'PropertyName1','Value1',...)
%    In this case, the assigment allows to link to internal or external links, 
%    as well as evaluated expression, with the syntax cases for the 'Value':
%     'field'                           a simple link to an other property 'field'
%     'field1.field2...'                a nested link to an other property
%     'file://some_file_path'           a local file URL
%     'http://some_distant_resource'    an HTTP URL (proxy settings may have to be set)
%     'https://some_distant_resource'   an HTTPS URL (proxy settings may have to be set)
%     'ftp://some_distant_resource'     an FTP URL (proxy settings may have to be set)
%     'matlab: some_expression'         some code to evaluate. 'this' refers to the 
%                                       object itself e.g. 'matlab: this.Signal*2'
%
%    File and URL can refer to compressed resources (zip, gz, tar, Z) which are 
%    extracted on-the-fly. In case the URL/file resource contains 'sections' a search token
%    can be specified with syntax such as 'file://filename#token'.
% 
%    SET(S) displays all object properties.
%
% Example: s=estruct; set(s, 'a', 1); set(s, 'b.c','a','alias'); get(s, 'b.c') == 1 && strcmp(get(s, 'b.c','alias'),'a')
%
% See also: fieldnames, findfield, isfield, get, estruct/getalias, estruct/get

% NOTE: the rationale here is to implement all the logic in subsref and just call it.

  if ~isa(s, 'estruct')
    builtin('set', s, varargin{:});
    return
  end
  
  if nargin==1, s = fieldnames(s); return; end

  % handle array of struct
  if numel(s) > 1
    for index=1:numel(s)
      s(index) = set(s(index), varargin{:});
    end
    return
  end
  
  % check last argument as 'alias' ? will not follow links, but set aliases directly.
  if nargin >=4 && rem(numel(varargin),2) == 1 && ischar(varargin{end}) && any(strcmp(varargin{end}, {'link','alias'}))
       follow = false;  % we set aliases and do not link for get/subsref.
       varargin(end) = [];
  else follow = true;   % we travel through links
  end
  
  % now 's' is a single object. We handle name/value pairs
  for index=1:2:numel(varargin)
    if index == numel(varargin), break; end
    name = varargin{index};
    value= varargin{index+1};
    if ~ischar(name) && ~iscellstr(name)
      error([ mfilename ': SET works with name/value pairs. The ' num2str(index) '-th argument is of type ' class(name) ]);
    end
    name = cellstr(name);
    for n_index=1:numel(name)
      if follow
        s = subsasgn(s, struct('type','.', 'subs',name{n_index}), value);
      else
        s = subsasgn(s, struct('type','()','subs',name{n_index}), value);
      end
    end
  end
  
  % reset cache (as we have changed the object: fields, values, ...)
  s.Private.cache.findfield = [];

