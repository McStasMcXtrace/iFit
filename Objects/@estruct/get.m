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
%    The PropertyName can be given as a cellstr, in which case all given properties
%    are returned.
%
%    V = GET(S,'PropertyName1.PropertyName2')
%    The 'PropertyName' can be a full structure path, such as 'field1.field2' in
%    in which case the value retrieval is made recursive.
%    When the retrieved value is itself a valid structure path (char), it is also 
%    travelled through, allowing 'aliases' such as in the following example:
%      s=estruct; set(s, 'a', 1); set(s, 'b.c','a'); 
%      get(s, 'b.c')  % returns s.a=1
%    This is equivalent to S('PropertyName') which follows links and PropertyName
%    can be a compound/complex property name.
%      get(s, 'b.c','nolink')  % returns s.b.c='a'
%
%    V = GET(S,'PropertyName','nolink') does not follows internal links
%    When the PropertyName points to a string value, it is returned as is.
% 
%    GET(S) displays all structure field names.
%
% Example: s=estruct; set(s, 'a', 1); set(s, 'b.c','a'); get(s, 'b.c') == 1 && strcmp(get(s, 'b.c','nolink'),'a')
%
% See also: fieldnames, findfield, isfield, set

  if ~isa(s, 'estruct')
    v = builtin('get', s, varargin{:});
    return
  end

  if nargin == 1, field=''; else field = varargin{1}; end
  if isempty(field), v = fieldnames(s); return; end
  if nargin <= 2, follow=true; else follow=false; end
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
  % get_single get a single field, recursively, and can follow char links

  if nargin < 3, follow=true; end
  % cut the field into pieces with '.' as separator
  [tok, rem] = strtok(field, '.');
  
  % get the highest level
  % subsref is faster than getfield which itself calls subsref
  v = subsref(s, struct('type','.','subs', tok)); 
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
    elseif strcmp(tok, 'matlab') && numel(v) > 8 % URL matlab:
      v = get_single_eval(s, v);
    elseif any(strcmp(tok, {'http' 'https' 'ftp' 'file'})) % URL http https ftp file: 
      try
        v = iLoad(v);
      end
    end
  end
  
  % ----------------------------------------------------------------------------
  function value = get_single_eval(this, value)
   % get_single_eval a sandbox to valuate a matlab expression
    try
      value = eval(value(8:end));
    catch ME
      disp(getReport(ME))
    end
  
