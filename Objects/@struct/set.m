function s = set(s, varargin)
% SET    Set structure properties.
%    V = SET(S,'PropertyName','Value') set the value of the specified
%    property/field in the structure.  
%    The 'PropertyName' can be a full structure path, such as 'field1.field2' in
%    in which case the value assigment is made recursive.
%
%    As opposed to GET(S) for structures, the assigment does not travel through  
%    valid aliases when the final value is a char.
% 
%    SET(S) displays all structure field names.

  if ~isstruct(s)
    builtin('set', s, varargin{:});
    return
  end
  field='';
  value=[];
  if nargin >=2,  field=varargin{1}; end
  if nargin >=3,  value=varargin{2}; end
  if isempty(field), s = fieldnames(s); return; end
  
  % handle array of struct
  if numel(s) > 1
    for index=1:numel(s)
      s(index) = set(s(index), field, value);
    end
    return
  end
  
  if ischar(field) && size(field, 1) > 1
    field = cellstr(field);
  end
  
  % handle array/cell of fields
  if iscellstr(field)
    for index=1:numel(field)
      s(index) = set(s, field{index}, value);
    end
    return
  end
  
  if ~ischar(field)
    error([ mfilename ': field to set in struct must be a char or cellstr, nor ' class(field) ]);
  end
  
  % cut the field into pieces with '.' as separator
  [tok, rem] = strtok(field, '.');
  
  if ~isfield(s, tok)
    s.(tok) = [];
  end
  
  % when rem is empty, we are were to set the value
  if isempty(rem)
    s = setfield(s, tok, value);
    return
  end
  
  % else get the sub-struct
  s2 = getfield(s, tok);
  
  % access deeper content recursively
  if ~isstruct(s2)
    s2 = []; % overwrite existing value
  end
  s2 = set(s2, rem(2:end), value);
  
  s = setfield(s, tok, s2); % update in parent struct
  
  if nargout == 0 && ~isempty(inputname(1)) && isa(s,'struct')
    assignin('caller',inputname(1),s);
  end
