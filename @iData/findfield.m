function [match, types, dims] = findfield(s, field, option)
% [match, types, nelements]=findfield(iData, field, option) : look for iData fields
%
%   @iData/findfield function to look for iData fields, type and number of elements
%
%   [match,type,n] = findfield(iData) returns the names of all iData fields
%   [match,type,n] = findfield(iData, field) returns the names of all iData fields 
%     that match 'field'
%   The optional 'option' may contain 'exact' to search for the exact occurence, and 'case'
%   to specifiy a case sensitive search.
%
% input:  s: object or array (iData)
%         field: field name to search, or '' (char).
%         option: 'exact' 'case' or '' (char)
% output: match: names of iData fields (cellstr)
%         types:  types of iData fields (cellstr), e.g. 'double', 'char', 'struct'...
%         nelements: total number of elements in iData fields (double)
% ex:     findfield(iData) or findfield(iData,'Title') or findfield(s,'Title','exact case')
%
% Version: $Revision: 1.8 $
% See also iData, iData/set, iData/get, iData/findobj, iData/findstr

% EF 23/09/07 iData implementation
% private function iData_getfields
% used in findstr and iData (find biggest field as the Signal)

if nargin == 1
  field = '';
end
if nargin <= 2
  option='';
end

if length(s(:)) > 1
  match = cell(1, length(s(:))); types=match; dims=match;
  for index=1:length(s)
    [m,t,n] = findfield(s(index), field, option);
    match{index}=m;
    types{index}=t;
    dims{index}=n;
  end
  return
end

[match, types, dims] = iData_getfields(struct(s), '');

if ~isempty(field)
  if isempty(strfind(option, 'case'))
    field= lower(field);
    matchs= lower(match);
  else
    matchs = match;
  end
  if strfind(option, 'exact')
    index = find(strcmp(field, matchs));
  else
    if iscellstr(field)
      index = [];
      for findex=1:length(field)
        tmp = strfind(matchs, field{findex});
        if iscell(tmp), tmp = find(cellfun('isempty', tmp) == 0); end
        index= [ index ; tmp ];
      end
      index = unique(index);
    else
      index = strfind(matchs, field);
    end
  end
  if ~isempty(index) && iscell(index), index = find(cellfun('isempty', index) == 0); end
  if isempty(index)
    match=[]; types=[]; dims=[];
  else
    match = match(index);
    types = types(index);
    dims  = dims(index);
  end
end

% ============================================================================
% private function iData_getfields, returns field, class, numel 
function [f, t, n] = iData_getfields(structure, parent)

f=[]; t=[]; n=[];
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
  n = cellfun('length', c);
end

if ~isempty(parent), f = strcat([ parent '.' ], f); end

% find sub-structures and make a recursive call for each of them
for index=transpose(find(cellfun('isclass', c, 'struct')))
  try
  [sf, st, sn] = iData_getfields(c{index}, f{index});
  f = [f(:) ; sf(:)];
  t = [t(:) ; st(:)];
  n = [n(:) ; sn(:)];
  end
end


