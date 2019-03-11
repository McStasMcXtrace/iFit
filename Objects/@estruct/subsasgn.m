function a = subsasgn(a,S,val)
% a = subsasgn(a,index,b) indexed assignement
%
%   @estruct/subsasgn: function defines indexed assignement 
%   such as a(1:2,3)    set indexed elements in array
%           a.('b')     set single field in structure, same as a.b
%           a('b')      set complex field in structure and follow links
%
%   The special syntax a('field') allows recursive field assigment and follows
%   links, e.g. a=estruct('f1',1,'f2','f1'); a('f2')=3; will set a.f1=3
%   and you may as well use: a('f3.a')=2 will set a.f3.a=1 directly.
%
% Version: $Date$

if numel(a) > 1
  a = builtin('subsasgn', a, S, val);
  return
end

% check for protected properties
if iscell(S(1).subs) && any(strcmp(S(1).subs{1}, a.Protected))
  error([ mfilename ': can not set Protected property ' S(1).subs{1} ' in object ' a.Tag ]);
elseif ischar(S(1).subs(1)) && any(strcmp(S(1).subs, a.Protected))
  error([ mfilename ': can not set Protected property ' S(1).subs ' in object ' a.Tag ]);
end

% we use the default subsasgn except for single level assigment
%   S.type = '()';
%   ischar(S.subs) -> set(a, S.subs, val)
if numel(S) == 1 && strcmp(S.type,'()') && iscellstr(S.subs) && ~strcmp(S.subs{1},':')
  set(a, S.subs{1}, val); % new syntax a('fields') = value supports links
  return
elseif numel(S) == 1 && strcmp(S.type,'.') && ischar(S.subs) && ~isfield(a, S.subs)
  set(a, S.subs, val, 'nolink');            % new field: a.('field') = value. does not follow links.
  return
end

% special syntax not handled: use builtin
a = builtin('subsasgn', a, S, val);

% reset cache (as we have changed the object: fields, values, ...)
a.Private.cache.findfield = [];
