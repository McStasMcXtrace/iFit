function [match, field] = findstr(s, str, option)
% FINDSTR look for strings stored in object
%   match = FINDSTR(s)      returns all strings in object
%
%   match = FINDSTR(s, str) returns the strings in object containg 'str'
%
%   match = FINDSTR(s, str, option) the 'option' may contain 'exact' to search 
%   for the exact occurence, and 'case' to specifiy a case sensitive search.
%
%   [match,field] = FINDSTR(..) also returns field name it appears in. 
%
% syntax:
%   [match, field]=findstr(s, str, option)
%
% input:  s:      object or array (struct)
%         str:    string to search in object, or '' (char or cellstr).
%         option: 'exact' 'case' or '' (char)
% output: match:  content of struct fields that contain 'str' (cellstr)
%         field:  name of struct fields that contain 'str' (cellstr)
%
% Example: s=estruct('x',1:10,'y','blah'); ischar(findstr(s, 'blah'))
% Version: $Date$ $
% See also estruct, set, get, findobj, findfield

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

[fields, types] = findfield(s,'','char');  % get all fields and types

index=[ find(strcmp('char', types)) ; find(strcmp('cell', types)) ]; 
if isempty(index), field=[]; match=[]; return; end

clear types
% get all field names containg char/cellstr data
fields = fields(index);
% get their content
matchs = get(s, fields);
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
