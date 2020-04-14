function a = subsasgn_single(a, S, val, a0)
% subsasgn_single single level assignment
%   SUBSASGN_SINGLE(a,S,val,a0) assigns a.(S) = val in root object a0
  if nargin<4, a0=a; end
  if ischar(S), S=struct('type','.','subs', S); end
  if ~isstruct(S),  error([ mfilename ': invalid reference (2nd arg) expect struct, is ' class(S) ]); end
  if numel(S) ~= 1, error([ mfilename ': only works with a single level reference' ]); end

  default = true;
  switch S.type
  case {'()','.'} % syntax: a('fields') does not follow links (setalias).
                  % syntax: a.('field') follows links (set), can be a compound field.
    if ischar(S.subs) S.subs = cellstr(S.subs); end
    if iscellstr(S.subs) && isscalar(S.subs)
      if ~isobject(a) || ~ismethod(a, S.subs{1})
        % follow links for '.' subsref, not for '()'
        a = set_single(a, S.subs{1}, val, S.type(1)=='.', a0);  % which handles aliases
        default = false;
      else  error([ mfilename ': ' class(a) '.' S.subs{1} ' is a reserved method name. Can not set new property.' ]);
      end
    end
  case '{}' % syntax: a{axis_rank} set axis value/alias (setaxis)
    if isa(a, 'estruct') && numel(S.subs{1}) == 1 % scalar numeric or char
      a = setaxis(a, S.subs{1}, val); % also set cache.check_requested to true
      default = false;
    end
  end
  if default  % other cases
    a = builtin('subsasgn', a, S, val);
  end

% ----------------------------------------------------------------------------
function s = set_single(s, field, value, follow, s0)
  % set_single set a single field to given value
  % when follow is true, the existing field value is checked for further link
  %   then the initial structure s0 is set again.

  if nargin <= 3, follow=true; end
  if nargin <= 4, s0=s; end

  if isa(s0, 'estruct') && s0.verbose > 2
    disp([ mfilename ': DEBUG: field ' field ' -> ' class(value) ' [' num2str(numel(value)) '] ' num2str(follow) ]);
  end

  % cut the field into pieces with '.' as separator
  if any(field == '.')
    field = textscan(field,'%s','Delimiter','.'); field=field{1};
    typs=cell(size(field)); [typs{:}] = deal('.');
    S = struct('type',typs, 'subs', field);
  else
    field = cellstr(field);
    S = struct('type','.','subs', field{1});
  end

  % use builtin subsasgn for the whole path when 'not follow'
  if ~follow && numel(S) > 1
    s = subsasgn(s, S, value);
    return
  else
    % now handle each level and test for char/alias value
    if numel(S) > 1 % oh no ! this is again recursive (multi levels)
      s = subsasgn_recursive(s, S, value, s0);
    else
      % single level indeed
      if ~isfield(s, S.subs) % new field ?
        if isa(s, 'estruct') && isempty(findprop(s, S.subs))
          s.addprop(S.subs);
        else
          s.(S.subs) = [];  % a normal structure
        end
        if isa(s0, 'estruct'), s0.Private.cache.findfield = []; end
      end
      % subsasgn is faster than setfield which itself calls subsasgn
      if follow
        % get possible alias
        v = builtin('subsref',s, S);
        if ischar(v) && isfield(s0, strtok(v, '.')) % link exists in original object ?
          try
            s0 = set_single(s0, v, value, follow); % set link in root object
            % must exit recursive levels
            return
          end
        elseif isa(s0, 'estruct') && isnumeric(v) && numel(v) ~= numel(value)
          s0.Private.cache.check_requested = true;
          s0.Private.cache.size = [];
          s0.Private.cache.std_w= [];
          s0.Private.cache.std_c= [];
        end
      end
      % if ~follow: change value. Calls "s.(tok)=value" i.e. subsasgn
      s = builtin('subsasgn',s, S, value); % set true value/alias (no follow)
    end
  end

