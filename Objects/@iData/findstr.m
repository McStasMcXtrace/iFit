function [match, field] = findstr(s, str, option)
% FINDSTR Look for strings within object.
%   match = FINDSTR(s)      returns all strings in object
%
%   match = FINDSTR(s, str) returns the strings in object containg 'str'
%
%   match = FINDSTR(s, str, option) the 'option' may contain 'exact' to search
%   for the exact occurence, 'case' to specifiy a case sensitive search,
%   'first' and 'last' to specify the first or last match.
%
%   [match,field] = FINDSTR(..) also returns field name it appears in.
%
% Example: s=iData('x',1:10,'y','blah'); ischar(findstr(s, 'blah'))
% Version: $Date$ $Version$ $Author$
% See also iData, set, get, findobj, findfield

if nargin == 1
  str='';
end
if nargin <= 2
  option='';
end

if numel(s) > 1
  match = cell(1, numel(s)); field=match;
  for index=1:numel(s)
    [m,f] = findstr(s(index), str, option);
    match{index}=m;
    field{index}=f;
  end
  return
end

% get all field names containing char/cellstr data
fields = findfield(s,'','char');

if isempty(fields), field=[]; match=[]; return; end

% get their content, but does not follow links
matchs = get(s, fields,'alias');
% convert to char using private inline function handle
matchs = cellfun(@cell2char, matchs, 'UniformOutput', false);

% all fields (no pattern given)
if isempty(str)
  match = matchs;
  field = fields;
  return
end

% multiple search
if iscellstr(str)
  match = cell(size(str)); field=match;
  for index=1:numel(str)
    [m,f] = struct_findstr_single(str{index}, option, fields, matchs);
    match{index}=m;
    field{index}=f;
  end
  return
end

% single search
[match, field] = struct_findstr_single(char(str), option, fields, matchs);

% -------------------------------------------------------------------------
function [match, field] = struct_findstr_single(str, option, field, match)

% search for the string
if strfind(option, 'exact')
  % handle 'exact' option
  index = strcmp(match, str);
elseif strfind(option, 'case')
  % handle 'case' option
  index = ~cellfun(@isempty, strfind(match, str));
else
  % relaxed search: non case sensitive, find token (not exact comparison)
  index = ~cellfun(@isempty, strfind(lower(match), lower(str)));
end
if ~isempty(strfind(option, 'first')) && ~isempty(index)
  index=index(1);
elseif ~isempty(strfind(option, 'last')) && ~isempty(index)
  index=index(end);
end
match = match(index);
field = field(index);

if numel(field) == 1, field=field{1}; end
if numel(match) == 1, match=match{1}; end

% -------------------------------------------------------------------------
function c=cell2char(c)
if iscell(c) && ischar(c{1})
    c = char(c);
elseif ~ischar(c)
    c = '';
end
c=c(:)';
