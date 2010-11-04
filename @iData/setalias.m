function s_out = setalias(a_in,names,links,labels)
% [s,...] = setalias(s, AliasName, AliasLink, AliasLabel) : set iData aliases
%
%   @iData/setalias function to set iData aliases.
%   The function works also when AliasName, AliasLink, AliasLabel
%     are given as cell strings. The AliasLink may be of any class, but char is
%     interpreted as a link to search in the object.
%   The special name 'this' may be used in Aliases to refer the object itself.
%   When the link is empty, the alias is removed, so that
%     setalias(s, alias)       deletes an alias
%     setalias(s, getalias(s)) deletes all alias definitions.
%   The command setalias(iData,'Signal') sets the Signal to the biggest numerical field.
%   The input iData object is updated if no output argument is specified.
%   An new field/alias may be defined with the quick syntax 's.alias = value'.
%
% input:  s: object or array (iData)
%         AliasName: Name of existing or new alias (char/cellstr)
%         AliasLink: definition of the alias, or '' to remove the alias (cell of char/double/...)
%         AliasLabel: optional description/label of the alias (char/cellstr)
% output: s: array (iData)
% ex:     setalias(iData,'Temperature','Data.Temperature','This is the temperature')
%         setalias(iData,'Temperature','this.Data.Temperature')
%         setalias(iData,'Temperature',1:20)
%         setalias(iData,'T_pi','[ this.Data.Temperature pi ]')
%
% Version: $Revision: 1.10 $
% See also iData, iData/getalias, iData/get, iData/set, iData/rmalias

% EF 27/07/00 creation
% EF 23/09/07 iData implementation

s_out=a_in;
if nargin == 1
  % makes a check of aliases, Signal, Error, Monitor, warns invalid ones.
  for index = 1:length(a_in(:))
    a = a_in(index); % current object in array/single element
    for j1=1:length(a.Alias.Names)
      try
        b=get(a, a.Alias.Names{j1});
      catch
        v = a.Alias.Values{j1};
        v = mat2str(v(1:min(20,length(v))));
        iData_private_warning(mfilename,[ 'the Alias ' a.Alias.Names{j1} '=' v ' is not valid in object ' inputname(1) ' ' a.Tag '.' ]);
      end
    end
  end
  return
elseif nargin == 2
  names = cellstr(names);
  if length(names) == 1 & strmatch(names{1}, 'Signal','exact')
    a_in = iData(a_in);
    return
  end
  links = ''; labels=''; % removes aliases
elseif nargin == 3
  labels='';
end

names = cellstr(names);
if ischar(links), links = cellstr(links); end
if ~iscell(links),links = { links }; end

labels= cellstr(labels);

s_out = a_in(:);

if isempty(names), return; end
for index = 1:length(s_out)
  a = s_out(index); % current object in array/single element
  
  for j1=1:length(names) % loop on alias names
    name = names{j1};
    if length(links) == length(names), link = links{j1}; else link=''; end
    if length(labels)== length(names), label= labels{j1}; else label=''; end
    
    % check that name is not a class member
    f = fieldnames(a);
    if ~isempty(strmatch(lower(name), lower(f), 'exact'))
      iData_private_error(mfilename,[ 'the Alias ' name ' is a protected name in object ' inputname(1) ' ' a.Tag '.' ]);
    end
    if ~isvarname(lower(name)) & isempty(findfield(a, name))
      iData_private_warning(mfilename,[ 'the Alias "' name '" is not a valid Alias name in object ' inputname(1) ' ' a.Tag '.' ]);
      continue
    end
    alias_names = a.Alias.Names; % this is a cellstr of Alias names
    alias_num   = strmatch(lower(name), lower(alias_names), 'exact');
    
    if isempty(link) & any(alias_num <= 3) 
      % set Signal, Error, Monitor to empty (default)
      a.Alias.Values{alias_num} = [];
      a.Alias.Labels{alias_num} = [];
    elseif isempty(link) & any(alias_num > 3) % protect Signal, Error, Monitor
      % remove these aliases from Alias list
      a.Alias.Names(alias_num)  = [];
      a.Alias.Values(alias_num) = [];
      a.Alias.Labels(alias_num) = [];
    elseif ~isempty(link)
      % update or add alias
      if isempty(alias_num) % add
        a.Alias.Names{end+1} = name;
        a.Alias.Values{end+1}= link;
        a.Alias.Labels{end+1}= regexprep(label,'\s+',' ');
      else
        a.Alias.Names{alias_num} = name;
        if ~isempty(link) | isempty(label), a.Alias.Values{alias_num}= link; end
        if ~isempty(label)| isempty(a.Alias.Labels{alias_num}), a.Alias.Labels{alias_num}= label; end
      end
    end
  end % for alias names
  a = iData_private_history(a, mfilename, a, name, link, label);
  
  s_out(index) = a;
end % for index

if length(s_out) > 1
  s_out = reshape(s_out,size(a_in));
end

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),s_out);
end
