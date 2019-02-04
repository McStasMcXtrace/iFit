function s = get(s, field)
% GET    Get structure properties.
%    V = GET(S,'PropertyName') returns the value of the specified
%    property/field in the structure.  If S is an array, then get will 
%    return a cell array of values.  
%    If 'PropertyName' is replaced by a 1-by-N or N-by-1 cell array of strings 
%    containing property names, then GET will return an M-by-N cell array of
%    values.
% 
%    GET(S) displays all structure field names.

  if nargin == 1, field=''; end
  if isempty(field), s = fieldnames(s); return; end
  
  % handle array of struct
  if numel(s) > 1
    sout = {};
    for index=1:numel(s)
      sout{end+1} = get(s(index), field);
    end
    sout = reshape(sout, size(s));
    s = sout; 
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
    s = sout; 
    return
  end
  
  if ~ischar(field)
    error([ mfilename ': field to search in struct must be a char or cellstr, nor ' class(field) ]);
  end
  
  % cut the field into pieces with '.' as separator
  [tok, rem] = strtok(field, '.');
  
  % get the highest level
  s = getfield(s, tok);
  if isempty(rem), return; end
  
  % an access deeper content recursively
  s = get(s, rem(2:end));
