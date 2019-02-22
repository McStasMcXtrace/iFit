function v = get(s, varargin)
% GET    Get structure properties.
%    V = GET(S,'PropertyName') returns the value of the specified
%    property/field in the structure.  If S is an array, then get will 
%    return a cell array of values.  
%    If 'PropertyName' is replaced by a 1-by-N or N-by-1 cell array of strings 
%    containing property names, then GET will return an M-by-N cell array of
%    values.
%    This is equivalent to S.('PropertyName') which does not follow links and
%    PropertyName must be a single property name.
%
%    V = GET(S,'PropertyName1.PropertyName2')
%    The 'PropertyName' can be a full structure path, such as 'field1.field2' in
%    in which case the value retrieval is made recursive.
%
%    V = GET(S,'PropertyName','link') follows internal links
%    When the retrieved value is itself a valid structure path (char), it is also 
%    travelled through, allowing 'aliases' such as in the following example:
%      s=estruct; set(s, 'a', 1); set(s, 'b.c','a'); 
%      get(s, 'b.c','link')  % returns s.a=1
%    This is equivalent to S('PropertyName') which follows links and PropertyName
%    can be a compound/complex property name.
% 
%    GET(S) displays all structure field names.
%
% See also: fieldnames, findfield, isfield, set

  if ~isa(s, 'estruct')
    v = builtin('get', s, varargin{:});
    return
  end

  if nargin == 1, field=''; else field = varargin{1}; end
  if isempty(field), v = fieldnames(s); return; end
  if nargin <= 2, follow=false; else follow=true; end
  v = [];
  
  if ischar(field) && size(field, 1) > 1
    field = cellstr(field);
  end
  
  % handle array of struct
  if numel(s) > 1
    sout = {};
    for index=1:numel(s)
      sout{end+1} = get(s(index), field, follow);
    end
    sout = reshape(sout, size(s));
    v = sout; 
    return
  end
  
  % handle array/cell of fields
  if iscellstr(field)
    sout = {};
    for index=1:numel(field)
      sout{end+1} = get(s, field{index}, follow);
    end
    sout = reshape(sout, size(field));
    v = sout; 
    return
  end
  
  if ~ischar(field)
    error([ mfilename ': field to search in struct must be a char or cellstr, nor ' class(field) ]);
  end
  
  v = get_single(s, field, follow);
  
% ------------------------------------------------------------------------------
function v = get_single(s, field, follow)

  if nargin < 3, follow=false; end
  % cut the field into pieces with '.' as separator
  [tok, rem] = strtok(field, '.');
  
  % get the highest level
  v = getfield(s, tok);
  if ~isempty(rem) && isstruct(v)
    % access deeper content recursively
    v = get_single(v, rem(2:end), follow);
  end
  
  % check if content is an alias, e.g. refers to a field/path
  if ischar(v) && follow
    % valid characters are A-Za-z_0-9
    tok = strtok(v, ' .()[]{};:=+-!@#$%^&*''"\|<>,?`~');
    if isfield(s, tok)
      try
        v = get_single(s, v, follow);
      end
    end
  end
  
