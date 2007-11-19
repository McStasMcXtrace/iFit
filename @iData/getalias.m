function [link, label, names] = getalias(s,alias)
% [link, label] = getalias(s, 'AliasName') : get iData alias
%
%   @iData/getalias function to get iData alias.
%   [link, label]          = getalias(s, alias) returns the alias link and its label/description.
%   [names, links, labels] = getalias(s)        returns all defined aliases.
%   The Signal, Error and Monitor aliases are always defined.
%
% input:  s: object or array (iData)
%         alias: alias name to inquire in object, or '' (char).
% output: link: alias link/definition (char/cellstr)
%         label: alias description (char/cellstr)
%         names: all defined alias names (cellstr)
% ex:     getalias(iData) or getalias(iData,'Signal')
%
% See also iData, iData/set, iData/get, iData/setalias

% EF 23/09/07 iData implementation
% ============================================================================

if nargin == 1
  alias = '';
end

if length(s(:)) > 1
  link = cell(size(s)); label=link; names=link;
  for index=1:length(s(:))
    [l,b,n] = getalias(s(index), alias);
    link{index} =l;
    label{index}=b;
    names{index}=n;
  end
  return
end

if isempty(alias)
  % NOTE: output arguments is shifted w.r.t. getalias(s,alias)
  link = s.Alias.Names(:);
  label= s.Alias.Values(:);
  names= s.Alias.Labels(:);
  return
end
names=alias;
alias_names = s.Alias.Names; % this is a cellstr of Alias names
alias_num   = strmatch(lower(alias), lower(alias_names), 'exact');

if isempty(alias_num)
  link=[]; label=[]; 
else
  link = s.Alias.Values{alias_num};
  label= s.Alias.Labels{alias_num};
end

