function s = set(s, varargin)
% SET    Set object properties.
%  V = SET(S,'PropertyName','Value') 
%    Set the value of the specified property/field in the structure.
%    The object S can be an array. The Value can be anything (char, alias, numeric).
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
%     'field'                         a simple link to an other property 'field'
%     'field1.field2...'              a nested link to an other property
%     'file://some_file_path'         a local file URL
%     'http://some_distant_resource'  an HTTP URL (proxy settings may need to be set)
%     'https://some_distant_resource' an HTTPS URL (proxy settings may need to be set)
%     'ftp://some_distant_resource'   an FTP URL (proxy settings may need to be set)
%     'matlab: some_expression'       some code to evaluate. 'this' refers to the 
%                                     object itself e.g. 'matlab: this.Signal*2'
%
%    File and URL can refer to compressed resources (zip, gz, tar, Z) which are 
%    extracted on-the-fly. In case the URL/file resource contains 'sections' a search token
%    can be specified with syntax such as 'file://filename#token'.
%
%  V = SET(S, 'Property','Value','Label')
%  Sets the Value and Label for the given Property.
% 
%  SET(S) displays all object properties.
%
%  To remove an alias/property use RMFIELD or RMALIAS.
%
% Example: s=iData; set(s, 'a', 1); set(s, 'b.c','a','alias'); get(s, 'b.c') == 1 && strcmp(get(s, 'b.c','alias'),'a')
% Version: $Date$ $Version$ $Author$
% See also iData, fieldnames, findfield, isfield, set, get, getalias, setalias, 
%   getaxis, setaxis, rmalias, rmfield

% NOTE: the rationale here is to implement all the logic in subasgn and just call it.

  if ~isa(s, 'iData')
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
  lab = ''; follow = true; % default: we travel through links
  % check last argument as 'alias' ? will not follow links, but set aliases directly.
  if nargin >=4 && rem(numel(varargin),2) == 1 && ischar(varargin{end}) 
    if any(strcmp(varargin{end}, {'link','alias'}))
       follow = false;  % we set aliases and do not link for get/subsref
    else lab = varargin(end);
    end
    varargin(end) = []; 
  end
  
  % now 's' is a single object. We handle name/value pairs
  for index=1:2:numel(varargin)
    if index == numel(varargin)
      warning([ mfilename ': SET works with name/value pairs. The ' num2str(index) '-th name argument misses its value. Skipping.' ]);
      break;
    end % missing value ?
    name = varargin{index};
    value= varargin{index+1};
    if ~ischar(name) && ~iscellstr(name) && ~(isscalar(name) && isnumeric(name))
      error([ mfilename ': SET works with name/value pairs. The ' num2str(index) '-th argument is of type ' class(name) ]);
    end
    if ischar(name), name = cellstr(name); end
    for n_index=1:numel(name)
      if isnumeric(name(n_index)) || ~isnan(str2double(name(n_index)))
        S.type = '{}'; S.subs={ name(n_index) };
        s = subsasgn(s, S, value);
      elseif follow
        s = subsasgn(s, struct('type','.', 'subs',name{n_index}), value);
      else
        s = subsasgn(s, struct('type','()','subs',name{n_index}), value);
      end
    end
  end
  
  if ~isempty(lab)
    label(s, name, lab);
  end
  
  % reset cache (as we have changed the object: fields, values, ...)
  s.Private.cache.findfield = [];

