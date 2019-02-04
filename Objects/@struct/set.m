function s = set(s, field, value)
% SET    Set structure properties.
%    V = SET(S,'PropertyName','Value') set the value of the specified
%    property/field in the structure.  
% 
%    SET(S) displays all structure field names.

  if nargin == 1, field=''; end
  if nargin < 2, value=[]; end
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
  
  % when rem is empty, we are were to set the value
  if isempty(rem)
    s = setfield(s, tok, value);
    return
  end
  
  % else get the sub-struct
  if ~isfield(s, tok)
    s.(tok) = [];
  end
  s2 = getfield(s, tok);
  
  % access deeper content recursively
  if ~isstruct(s2)
    s2 = []; % overwrite existing value
    s2.(rem(2:end)) = value;
  else
    s2 = set(s2, rem(2:end), value);
  end
  
  s = setfield(s, tok, s2); % update in parent struct
