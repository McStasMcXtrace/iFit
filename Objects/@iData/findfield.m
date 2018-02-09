function [match, types, dims] = findfield(s, field, option)
% [match, types, nelements]=findfield(iData, field, option) : look for iData fields
%
%   @iData/findfield function to look for iData fields, type and number of elements
%
%   [match,type,n] = findfield(iData) returns the names of all iData fields
%   [match,type,n] = findfield(iData, field) returns the names of all iData fields 
%     that match 'field'
%   The optional 'option' argument can contain one or more keywords:
%     The 'exact'   option will search for the exact occurences.
%     The 'case'    option specifies a case sensitive search. 
%     The 'numeric' option will return only numerical fields.
%     The 'char'    option will return only character fields.
%     The 'cache'   option allows to re-use a previous field search cache for  
%                     faster consecutive searches on the same object.
%     The 'first'   option returns the first match (with shortest name/link).
%     The 'biggest' option returns the biggest match (has max nb of elements).
%
%   The cache mechanism allows to do iterative searches on the same object, e.g.
%     findfield(s,'Title') and then findfield(s,'Signal','cache numeric')
%
% input:  s: object or array (iData)
%         field: field name to search, or '' (char).
%         option: empty or 'exact' 'case' 'char' or 'numeric' (char)
% output: match:  names of iData fields (cellstr)
%         types:  types of iData fields (cellstr), e.g. 'double', 'char', 'struct'...
%         nelements: total number of elements in iData fields (double)
% ex:     findfield(iData) or findfield(iData,'Title') or findfield(s,'Title','exact case')
%
% Version: $Date$
% See also iData, iData/set, iData/get, iData/findobj, iData/findstr

% EF 23/09/07 iData implementation
% private function iData_getfields
% used in findstr and iData (find biggest field as the Signal)

persistent cache

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

% test if the cache can be used: requires Tag to be the same
if ~isempty(cache) && ~strcmp(cache.Tag , s.Tag)
  cache = [];
end

if isempty(cache) || isempty(strfind(option, 'cache'))

  warning off MATLAB:structOnObject
  
  struct_s=struct(s);
  struct_s=rmfield(struct_s,'Alias');
  struct_s=rmfield(struct_s,'Command');
  struct_s=rmfield(struct_s,'Tag');
  % add the Aliases
  aliases=getalias(s);

  for index=1:numel(aliases)
      alias=aliases{index};
      value = get(s,alias);
      if isa(value, 'iData') || isa(value, 'iFunc')
        value = struct(value);
      end
      struct_s.(alias) = value;
  end

  % find all fields in object structure
  [match, types, dims] = iData_getfields(struct_s, '');

  % add fields from Alias content (when numeric/structure)
  for index=getalias(s)'
    this = getalias(s,index{1}); % alias def or value
    if isstruct(this)
      [sf, st, sn] = iData_getfields(this, index{1});
      match = [match(:) ; sf(:)];
      types = [types(:) ; st(:)];
      dims  = [dims(:) ;  sn(:)];
      clear sf st sn
    elseif isnumeric(this)
      match{end+1} = index{1};
      types{end+1} = class(this);
      dims(end+1)  = numel(this);
    end
  end
  [match, index] = unique(match);
  types=types(index);
  dims =dims(index);
  % store [match, types, dims] in cache
  cache.match = match;
  cache.types = types;
  cache.dims  = dims;
  cache.time  = clock;
  cache.Tag   = s.Tag;
else % restore from cache
  match = cache.match;
  types = cache.types;
  dims  = cache.dims;
end


% filter fields: numeric and char types
if ~isempty(strfind(option, 'numeric')) || ~isempty(strfind(option, 'char'))
  % remove fields that we do not want as Signal
  index=[ find(strcmp('Date', match)) find(strcmp('ModificationDate', match)) ] ;
  types(index) = {'char'};
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
    for findex=1:length(field)
      % get direct name or from last word of 'matchs'
      m     = find_last_word(matchs);
      this_index = [ find(strcmp(field{findex}, matchs)) ; find(strcmp(field{findex}, m)) ];
      index = [ index this_index ];
    end
    index = unique(index);
  else
    if iscell(field) && ischar(field{1})
      index = [];
      for findex=1:length(field)
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
% private function iData_getfields, returns field, class, numel 
function [f, t, n] = iData_getfields(structure, parent)

f={}; t={}; n=[];
if ~isstruct(structure), return; end

if numel(structure) > 1
  structure=structure(:);
  for index=1:length(structure)
    [sf, st, sn] = iData_getfields(structure(index), [ parent '(' num2str(index) ')' ]);
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
  [sf, st, sn] = iData_getfields(c{index}, f{index});
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
