function [match, field] = findstr(s, str, option)
% [match, field]=findstr(s, str, option) : look for strings stored in iData
%
%   @iData/findstr function to look for strings stored in iData
%
%   [match,field] = findfstr(iData, field) returns the string containg 'str' 
%     and the field name it appears in. If 'str' is set to '', the content of all
%     character fields is returned.
%   The 'option' may contain 'exact' to search for the exact occurence, and 'case'
%   to specifiy a case sensitive search.
%
% input:  s: object or array (iData)
%         str: string to search in object, or '' (char).
%         option: 'exact' 'case' or '' (char)
% output: match: content of iData fields that contain 'str' (cellstr)
%         field: name of iData fields that contain 'str' (cellstr)
% ex:     findstr(iData,'ILL') or findstr(s,'TITLE','exact case')
%
% Version: $Revision: 1.8 $
% See also iData, iData/set, iData/get, iData/findobj, iData/findfield

% EF 23/09/07 iData implementation

if nargin == 1
  str='';
end
if nargin <= 2
  option='';
end

if length(s(:)) > 1
  s=s(:);
  match = cell(1, length(s)); field=match;
  for index=1:length(s)
    [m,f] = findstr(s(index), str, option);
    match{index}=m;
    field{index}=f;
  end
  return
end

field={}; match={};
[fields, types] = findfield(s);  % get all fields and types

index=find(strcmp('char', types));
if isempty(index), match=[]; return; end

fields = fields(index); % get all field names containg char data

for index=1:length(fields)
  strfield = get(s, fields{index}); % get char content of field
  if ~ischar(strfield), continue; end
  if isempty(strfind(option, 'case'))
    str = lower(str);
    strs= lower(strfield);
  else
    strs= strfield;
  end
  % search for the pattern
  if isempty(str), j1=1;
  else
    if strfind(option, 'exact')
      j1 = find(strcmp(str, strs));
    else
      j1 = [ find(strncmp(str, strs,length(str))) strfind(strs, str) ]; 
    end
  end
  if ~isempty(j1)
    match = [ match ; strfield ];
    field = [ field ; fields{index} ];
  end
end

