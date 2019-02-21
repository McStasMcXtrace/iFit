function v = get(s, varargin)
% GET    Get structure properties.
%    V = GET(S,'PropertyName') returns the value of the specified
%    property/field in the structure.  If S is an array, then get will 
%    return a cell array of values.  
%    If 'PropertyName' is replaced by a 1-by-N or N-by-1 cell array of strings 
%    containing property names, then GET will return an M-by-N cell array of
%    values.
%
%    The 'PropertyName' can be a full structure path, such as 'field1.field2' in
%    in which case the value retrieval is made recursive.
%    When the retrieved value is itself a valid structure path (char), it is also 
%    travelled through, allowing 'aliases' such as in the following example:
%      s=estruct; set(s, 'a', 1); set(s, 'b.c','a'); 
%      get(s, 'b.c')  % returns s.a=1
% 
%    GET(S) displays all structure field names.

  if ~isa(s, 'estruct')
    v = builtin('get', s, varargin{:});
    return
  end

  if nargin == 1, field=''; else field = varargin{1}; end
  if isempty(field), v = fieldnames(s); return; end
  v = [];
  
  % handle array of struct
  if numel(s) > 1
    sout = {};
    for index=1:numel(s)
      sout{end+1} = get(s(index), field);
    end
    sout = reshape(sout, size(s));
    v = sout; 
    return
  end
  
  if ischar(field) && size(field, 1) > 1
    field = cellstr(field);
  end
  
  % handle array/cell of fields
  if iscellstr(field)
    sout = {};
    for index=1:numel(field)
      sout{end+1} = get(s, field{index});
    end
    sout = reshape(sout, size(field));
    v = sout; 
    return
  end
  
  if ~ischar(field)
    error([ mfilename ': field to search in struct must be a char or cellstr, nor ' class(field) ]);
  end
  
  v = get_single(s, field);
  
% ------------------------------------------------------------------------------
function v = get_single(s, field)

  % cut the field into pieces with '.' as separator
  [tok, rem] = strtok(field, '.');
  
  % get the highest level
  v = getfield(s, tok);
  if ~isempty(rem) && isstruct(v)
    % an access deeper content recursively
    v = get_single(v, rem(2:end));
  end
  
  % check if content is an alias, e.g. refers to a field/path
  if ischar(v)
    % valid characters are A-Za-z_0-9
    tok = strtok(v, ' .()[]{};:=+-!@#$%^&*''"\|<>,?`~');
    if isfield(s, tok)
      try
        v = get_single(s, v);
      end
    end
  end
  
