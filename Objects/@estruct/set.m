function s = set(s, varargin)
% SET    Set structure properties.
%    V = SET(S,'PropertyName','Value') set the value of the specified
%    property/field in the structure.  
%
%    V = SET(S,'PropertyName1.PropertyName2', 'Value')
%    The 'PropertyName' can be a full structure path, such as 'field1.field2' in
%    in which case the value assigment is made recursive.
%    When the target property is itself a valid structure path (char), it is also 
%    travelled through before assigment.
%
%    V = SET(S,'PropertyName1.PropertyName2', 'Value','nolink')
%    When the PropertyName points to a string value, it is assigned without 
%    travelling through it.
% 
%    SET(S) displays all structure field names.
%
% Example: s=estruct; set(s, 'a', 1); set(s, 'b.c','a'); get(s, 'b.c') == 1 && strcmp(get(s, 'b.c','nolink'),'a')
%
% See also: fieldnames, findfield, isfield, get

  if ~isa(s, 'estruct')
    builtin('set', s, varargin{:});
    return
  end
  
  field=''; value=[]; follow=false;
  if nargin >=2,  field=varargin{1}; end
  if isempty(field), s = fieldnames(s); return; end
  if nargin >=3,  value=varargin{2}; end
  if nargin >=4,  follow=true; end
  
  if ischar(field) && size(field, 1) > 1
    field = cellstr(field);
  end
  
  % handle array of struct
  if numel(s) > 1
    for index=1:numel(s)
      s(index) = set(s(index), field, value, follow);
    end
    return
  end
  
  % handle array/cell of fields
  if iscellstr(field)
    for index=1:numel(field)
      s(index) = set(s, field{index}, value, follow);
    end
    return
  end
  
  if ~ischar(field)
    error([ mfilename ': field to set in struct must be a char or cellstr, nor ' class(field) ]);
  end
  
  s = set_single(s, field, value, follow, s);
  
  % reset cache (as we have changed the object: fields, values, ...)
  s.Private.cache.findfield = [];
  
% ----------------------------------------------------------------------------
function [s, rec] = set_single(s, field, value, follow, s0)
  % set_single set a single field to given value
  % when follow is true, the existing field value is checked for further link
  %   then the initial structure s0 is set again.
  %
  % when the returned argument rec is true, recursivity is short-cut.
  if nargin < 5, s0=s; end
  if nargin < 4, follow=true; end
  rec=true;
  % cut the field into pieces with '.' as separator
  [tok, rem] = strtok(field, '.');
  
  if ~isfield(s, tok) % new field ?
    if isa(s, 'estruct')
      s.addprop(genvarname(tok));
    else
      s.(tok) = [];
    end
  end
  
  % when rem is empty, we are were to set the value
  if isempty(rem)
    if follow % check if the existing value is a char and valid link
      v = subsref(s, struct('type','.','subs', tok)); % do not follow, get possible alias
      target = strtok(v, ' .()[]{};:=+-!@#$%^&*''"\|<>,?`~');
      if isfield(s0, target)
        try
          s = set_single(s0, target, value, follow); % follow link
          rec=false;  % break recursivity to go back directly to root
          return
        end
      end
    end
    % change value. Calls "s.(tok)=value" i.e. subsasgn
    s = subsasgn(s, struct('type','.','subs', tok), value);
    return
  end
  
  % else get the sub-struct
  % subsref is faster than getfield which itself calls subsref
  s2 = subsref(s, struct('type','.','subs', tok));
  
  % access deeper content recursively
  if ~isstruct(s2)
    s2 = struct(); % overwrite existing value
  end
  
  [s2, rec] = set_single(s2, rem(2:end), value, follow, s0); % recursive
  if rec
    % update in parent struct. Calls "s.(tok)" i.e. subsref
    s = subsasgn(s, struct('type','.','subs', tok), s2);
  else
    s = s2; % we handled a link, and directly return the root object
  end
  
