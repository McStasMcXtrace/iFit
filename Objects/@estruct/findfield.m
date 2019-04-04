function [match, types, dims, sz] = findfield(s, field, option)
% FINDFIELD find fields, type and number of elements in object
%   match = FINDFIELD(s) returns the name of all fields
%
%   match = FINDFIELD(s, field)  returns the name of all fields that match 'field'
%
%   [match,type,n] = FINDFIELD(s, field, option)
%     The optional 'option' argument can contain one or more keywords:
%     The 'exact'   option will search for the exact occurrences.
%     The 'case'    option specifies a case sensitive search. 
%     The 'numeric' option will return only numerical fields.
%     The 'char'    option will return only character fields.
%     The 'first'   option returns the first match (with shortest name/link).
%     The 'biggest' option returns the biggest match (has max nb of elements).
%     The 'force'   option explicitly request a full scan of the object.
%     For instance, to identify the largest numeric element in object, use:
%       f=findfield(s,'','numeric biggest')
%
%   [match,type,n,sz] = FINDFIELD(...) also return the type, number of elements 
%   and size of fields
%
% syntax:
%   [match, types, nelements,sz]=findfield(s, field, option)
%
% input:  s: object or array (s)
%         field: field name to search, or '' (char).
%         option: empty or 'exact' 'case' 'char' or 'numeric' (char)
% output: match:  names of matching fields (cellstr)
%         types:  types of matching fields (cellstr), e.g. 'double', 'char', 'struct'...
%         nelements: total number of elements of matching fields (double)
%         sz:     size of matching fields (cellstr)
%
% Example: s=estruct('x', rand(6)); f=findfield(s,'','biggest'); numel(get(s,f)) == 36
% Version: $Date$ (c) E.Farhi. License: EUPL.
% See also estruct, estruct/set, estruct/get, estruct/findstr

if nargin == 1
  field = '';
end
if nargin <= 2
  option='';
end

if numel(s) > 1
  match = cell(1, numel(s)); types=match; dims=match;
  for index=1:numel(s)
    [m,t,n,z] = findfield(s(index), field, option);
    match{index}=m;
    types{index}=t;
    dims{index} =n;
    sz{index}   =z;
  end
  return
end



if ~isempty(strfind(option, 'force'))
  s.Private.cache.(mfilename) = [];
end

if isfield(s.Private, 'cache') && isfield(s.Private.cache, mfilename) ...
  && ~isempty(s.Private.cache.(mfilename))
  % get content from cache
  cache = s.Private.cache.(mfilename);
  match = cache.match;
  types = cache.types;
  dims  = cache.dims;
  sz    = cache.sizes;
else
  struct_s=rmfield(struct(s),s.Protected);  % private/protected stuff must remain so.
  
  % find all fields in object structure
  [match, types, dims, sz] = struct_getfields(struct_s, '');
  
  [match, index] = unique(match);
  types=types(index);
  dims =dims(index);
  sz   =sz(index);
  
  % set content into cache
  cache.match = match;
  cache.types = types;
  cache.dims  = dims;
  cache.sizes = sz;
  s.Private.cache.(mfilename) = cache;
end

% filter fields: numeric and char types
if ~isempty(strfind(option, 'numeric')) || ~isempty(strfind(option, 'char'))

  % now get the numeric ones
  if ~isempty(strfind(option, 'numeric'))   % e.g. double, int, ...
    index=          find(strcmp( 'double', types));
    index=[ index ; find(strcmp( 'single', types)) ];
    index=[ index ; find(strcmp( 'logical',types)) ];
    index=[ index ; find(strncmp('uint',   types, 4)) ];
    index=[ index ; find(strncmp('int',    types, 3)) ];
  elseif ~isempty(strfind(option, 'char'))  % 'char' option
    index=          find(strcmp( 'char', types));
  else index=[];                            % all / other
  end
  
  if ~isempty(index)
    match  = match(index); % restrict search to numeric/char data
    dims   = dims(index);
    types  = types(index);
    sz     = sz(index);
  end
  
  % sort fields in descending size order
  [dims, index]  = sort(dims, 'descend');
  match  = match(index);
  types  = types(index);
  sz     = sz(index);
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
  % get only the non-empty fields
  if ~isempty(index) && iscell(index), index = find(cellfun(@isempty, index) == 0); end
  if isempty(index)
    match={}; types={}; dims=[]; sz={};
  else
    match = match(index);
    types = types(index);
    dims  = dims(index);
    sz    = sz(index);
  end

  if isempty(match), return; end
end % field search on token

if ~isempty([ strfind(option, 'first') strfind(option, 'shortest') strfind(option, 'simplest') ])
  % we select the 'shortest' match
  [m, index] = min(cellfun(@length, match));
  match = match{index};
  types = types{index};
  dims  = dims(index);
  sz    = sz{index};
elseif ~isempty([ strfind(option, 'biggest') strfind(option, 'largest') ])
  % we select the 'biggest' match
  [m, index] = max(dims);
  match = match{index};
  types = types{index};
  dims  = dims(index);
  sz    = sz{index};
end


% ============================================================================
% private function struct_getfields, returns field, class, numel 
function [f, t, n, z] = struct_getfields(structure, parent)

f={}; t={}; n=[]; z={};
if ~isstruct(structure), return; end

% handle array of structures
if numel(structure) > 1
  structure=structure(:);
  for index=1:length(structure)
    [sf, st, sn, sz] = struct_getfields(structure(index), [ parent '(' num2str(index) ')' ]);
    f = [f(:) ; sf(:)];
    t = [t(:) ; st(:)];
    n = [n(:) ; sn(:)];
    z = [z(:) ; sz(:)];
  end
  return
end


f = fieldnames(structure);
f0= f;  % so that we keep initial single level fields

if ~isempty(parent) && ~isempty(f), % assemble compound field path
  % f = strcat([ parent '.' ], f); % this is time consuming
  for ii=1:length(f)
    f{ii} = [ parent '.' f{ii} ]; % slightly faster
  end
end

% get info by reading the whole content (when smaller than 1Gb, but fast)
failed=true;
if getfield(whos('structure'), 'bytes') < 1e9
  try
    c = struct2cell(structure); 
    t = cellfun(@class, c, 'UniformOutput', 0);
    n = cellfun(@numel, c);
    z = cellfun(@size,  c, 'UniformOutput', 0);
    failed = false; % get there: all is fine
    clear c
  end
end
% when storage is too large or failed, go through all fields sequentially (slower)
if failed
  t = cell(1, numel(f)); n=zeros(1,numel(f)); z=t;
  for ii=1:length(f)
    t{ii} = class(builtin('subsref',structure, struct('type','.','subs',f{ii})));
    z{ii} = size( builtin('subsref',structure, struct('type','.','subs',f{ii})));
    n(ii) = prod(z{ii});
  end
end

% find sub-structures and make a recursive call for each of them
for ii=1:length(f)
  if any(strcmp(t{ii},{'struct','estruct'}))  % structure -> recursive
    try
      [sf, st, sn,sz] = struct_getfields(...
        builtin('subsref',structure, struct('type','.','subs', f0{ii})), f{ii});
      f = [f(:) ; sf(:)];
      t = [t(:) ; st(:)];
      n = [n(:) ; sn(:)];
      z = [z(:) ; sz(:)];
    end
  end
end

% ------------------------------------------------------------------------------
function m = find_last_word(matchs)
% find_last_word get the last word in a cellstr

% m=strtrim(cellstr(fliplr(char(strtok(cellstr(fliplr(char(matchs))),'. ')))));
  m = cellfun(@(s)s((find(s == '.', 1, 'last')+1):end), matchs, 'UniformOutput', false);
  index = cellfun(@isempty, m); m(index) = matchs(index);
