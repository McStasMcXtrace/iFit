function [match, types, dims] = findfield(s, field, option)
% [match, types, nelements]=findfield(struct, field, option) : look for structure fields
%
%   @struct/findfield function to look for struct fields, type and number of elements
%
%   [match,type,n] = findfield(struct) returns the names of all struct fields
%   [match,type,n] = findfield(struct, field) returns the names of all struct fields 
%     that match 'field'
%   The optional 'option' argument can contain one or more keywords:
%     The 'exact'   option will search for the exact occurences.
%     The 'case'    option specifies a case sensitive search. 
%     The 'numeric' option will return only numerical fields.
%     The 'char'    option will return only character fields.
%     The 'first'   option returns the first match (with shortest name/link).
%     The 'biggest' option returns the biggest match (has max nb of elements).
%
% input:  s: object or array (struct)
%         field: field name to search, or '' (char).
%         option: empty or 'exact' 'case' 'char' or 'numeric' (char)
% output: match:  names of struct fields (cellstr)
%         types:  types of struct fields (cellstr), e.g. 'double', 'char', 'struct'...
%         nelements: total number of elements in struct fields (double)
% ex:     findfield(struct) or findfield(struct,'Title') or findfield(s,'Title','exact case')
%
% Version: $Date$ $Version$ $Author$
% See also struct, struct/set, struct/get, struct/findstr

if nargin == 1
  field = '';
end
if nargin <= 2
  option='';
end

if numel(s) > 1
  match = cell(1, numel(s)); types=match; dims=match;
  for index=1:numel(s)
    [m,t,n] = findfield(s(index), field, option);
    match{index}=m;
    types{index}=t;
    dims{index}=n;
  end
  return
end

struct_s=struct(s);

% find all fields in object structure
[match, types, dims] = struct_getfields(struct_s, '');

[match, index] = unique(match);
types=types(index);
dims =dims(index);

% filter fields: numeric and char types
if ~isempty(strfind(option, 'numeric')) || ~isempty(strfind(option, 'char'))

  % now get the numeric ones
  if ~isempty(strfind(option, 'numeric'))
    index=          find(strcmp( 'double', types));
    index=[ index ; find(strcmp( 'single', types)) ];
    index=[ index ; find(strcmp( 'logical',types)) ];
    index=[ index ; find(strncmp('uint',   types, 4)) ];
    index=[ index ; find(strncmp('int',    types, 3)) ];
  else % 'char' option
    index=          find(strcmp( 'char', types));
  end
  
  match  = match(index); % get all field names containing double data
  dims   = dims(index);
  types  = types(index);
  
  % sort fields in descending size order
  [dims, index]  = sort(dims, 'descend');
  match  = match(index);
  types  = types(index);
end

% restrict search to a given token 
if ~isempty(field)
  if isempty(strfind(option, 'case'))
    field = lower(field);
    matchs= lower(match);
  else
    matchs = match;
  end
  if strfind(option, 'exact')
    field = cellstr(field);
    index = [];
    for findex=1:numel(field)
      % get direct name or from last word of 'matchs'
      m     = find_last_word(matchs);
      this_index = [ find(strcmp(field{findex}, matchs)) ; find(strcmp(field{findex}, m)) ];
      index = [ index this_index ];
    end
    index = unique(index);
  else
    if iscell(field) && ischar(field{1})
      index = [];
      for findex=1:numel(field)
        tmp = strfind(matchs, field{findex});
        if iscell(tmp), tmp = find(cellfun('isempty', tmp) == 0); end
        index= [ index ; tmp ];
      end
      index = unique(index);
    else
      index = strfind(matchs, field);   % faster
    end
  end
  if ~isempty(index) && iscell(index), index = find(cellfun('isempty', index) == 0); end
  if isempty(index)
    match={}; types={}; dims=[];
  else
    match = match(index);
    types = types(index);
    dims  = dims(index);
  end

  if isempty(match), return; end
  if strfind(option, 'first')
    % we select the 'shortest' match
    [m, index] = min(cellfun(@length, match));
    match = match{index};
    types = types{index};
    dims  = dims(index);
  elseif strfind(option, 'biggest')
    % we select the 'biggest' match
    [m, index] = max(dims);
    match = match{index};
    types = types{index};
    dims  = dims(index);
  end
end

% ============================================================================
% private function struct_getfields, returns field, class, numel 
function [f, t, n] = struct_getfields(structure, parent)

f={}; t={}; n=[];
if ~isstruct(structure), return; end

if numel(structure) > 1
  structure=structure(:);
  for index=1:length(structure)
    [sf, st, sn] = struct_getfields(structure(index), [ parent '(' num2str(index) ')' ]);
    % sf = strcat(sf, [ '(' num2str(index) ')' ])
    f = [f(:) ; sf(:)];
    t = [t(:) ; st(:)];
    n = [n(:) ; sn(:)];
  end
  return
end

c = struct2cell(structure);
f = fieldnames(structure);
try
  t = cellfun(@class, c, 'UniformOutput', 0);
  n = cellfun(@numel, c);
catch
  t=cell(1,length(c));
  index = cellfun('isclass', c, 'double'); t(find(index)) = {'double'};
  index = cellfun('isclass', c, 'single'); t(find(index)) = {'single'};
  index = cellfun('isclass', c, 'logical');t(find(index)) = {'logical'};
  index = cellfun('isclass', c, 'char');   t(find(index)) = {'char'};
  index = cellfun('isclass', c, 'struct'); t(find(index)) = {'struct'};
  index = cellfun('isclass', c, 'cell');   t(find(index)) = {'cell'};
  index = cellfun('isclass', c, 'uint8');  t(find(index)) = {'uint8'};
  index = cellfun('isclass', c, 'uint16'); t(find(index)) = {'uint16'};
  index = cellfun('isclass', c, 'uint32'); t(find(index)) = {'uint32'};
  index = cellfun('isclass', c, 'uint64'); t(find(index)) = {'uint64'};
  index = cellfun('isclass', c, 'int8');   t(find(index)) = {'int8'};
  index = cellfun('isclass', c, 'int16');  t(find(index)) = {'int16'};
  index = cellfun('isclass', c, 'int32');  t(find(index)) = {'int32'};
  index = cellfun('isclass', c, 'int64');  t(find(index)) = {'int64'};
  n = cellfun('prodofsize', c);
end

if ~isempty(parent) && ~isempty(f), 
  % f = strcat([ parent '.' ], f); % this is time consuming
  for ii=1:length(f)
    f{ii} = [ parent '.' f{ii} ]; % slightly faster
  end
end

% find sub-structures and make a recursive call for each of them
for index=transpose(find(cellfun('isclass', c, 'struct')))
  try
  [sf, st, sn] = struct_getfields(c{index}, f{index});
  f = [f(:) ; sf(:)];
  t = [t(:) ; st(:)];
  n = [n(:) ; sn(:)];
  end
end

% ------------------------------------------------------------------------------
function m = find_last_word(matchs)
% m=strtrim(cellstr(fliplr(char(strtok(cellstr(fliplr(char(matchs))),'. ')))));
  m = cellfun(@(s)s((find(s == '.', 1, 'last')+1):end), matchs, 'UniformOutput', false);
  index = cellfun(@isempty, m); m(index) = matchs(index);
