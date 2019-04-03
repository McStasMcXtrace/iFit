function v = subsref_single(v, S, a)
% subsref_single single level subsref
  if nargin<3, a=v; end
  if ischar(S), S=struct('type','.','subs', S); end
  if ~isstruct(S),  error([ mfilename ': invalid reference (2nd arg) expect struct, is ' class(S) ]); end
  if numel(S) ~= 1, error([ mfilename ': only works with a single level reference' ]); end
  
  default = true;
  switch S.type
  case {'()','.'} % syntax: a('fields') does not follow links (getalias).
                  % syntax: a.('field') follows links (get), can be a compound field.
    if ischar(S.subs) S.subs = cellstr(S.subs); end
    if iscellstr(S.subs) && isscalar(S.subs) && ~ismethod(v, S.subs{1})
      % follow links for '.' subsref, not for '()'
      v = get_single(v, S.subs{1}, S.type(1)=='.', a);  % which handles aliases in 'a'
      default = false;
    end
  case '{}' % syntax: a{axis_rank} get axis value/alias (getaxis)
    if isa(v, 'estruct') && numel(S.subs{1}) == 1 % scalar numeric or char
      v = getaxis(v, S.subs{1});
      default = false;
    end
  end
  if default  % other cases
    v = builtin('subsref', v, S);
  end

% ------------------------------------------------------------------------------
function v = get_single(s, field, follow, s0)
  % get_single get a single field, recursively, and can follow char links
  % when follow is true, the existing field value is checked for further link
  %   then the initial structure s0 is set again.

  if nargin < 3, follow=true; end
  if nargin < 4, s0=s; end
  
  % cut the field into pieces with '.' as separator
  if any(field == '.')
    field = textscan(field,'%s','Delimiter','.'); field=field{1};
    typs=cell(size(field)); [typs{:}] = deal('.');
    S = struct('type',typs, 'subs', field);
  else
    field = cellstr(field);
    S = struct('type','.','subs', field{1});
  end
  
  % use builtin subsref for the whole path when 'not follow'
  if ~follow && numel(S) > 1
    v = builtin('subsref', s, S);
    return
  else
    % now handle each level and test for char/alias value
    v = s;
    for index=1:numel(S)
      % subsref is faster than getfield which itself calls subsref
      v = builtin('subsref',v, S(index)); % get true value/alias (no follow)
      % when 'follow', we need to access values iteratively and check for possible 'alias' (char)
      if follow && ischar(v) && ~strcmp(S(index).subs,'Source')
        if any(v == '.') && all(isletter(v) | v == '_' | v == '.' | isstrprop(v, 'digit')) % a compound link
          v = get_single(s, v, follow, s0);
        else
          v = get_single_alias(s0, v); % access a link/alias in initial object/structure
        end
      end
    end
  end
  
% ----------------------------------------------------------------------------
function v = get_single_alias(s, v)
  if ~ischar(v), return; end
  if isfield(s, v) % a link to a valid field in root object (single level)
    v = builtin('subsref',s, struct('type','.','subs', v)); % get true value/alias (no follow)
  elseif strcmp(v, 'matlab: sqrt(this.Signal)') % for default Error (no eval)
    v = sqrt(get(s,'Signal'));
  elseif strncmp(v, 'matlab',6) && numel(v) > 8 % URL matlab:
    v = get_single_eval(s, v);
  elseif ~isempty(dir(v)) || any(strcmp(strtok(v, ':'), {'http' 'https' 'ftp' 'file'})) % URL http https ftp file: 
    try
      v = iLoad(v);
    end
  end
  
% ----------------------------------------------------------------------------
function value = get_single_eval(this, value)
 % get_single_eval a sandbox to valuate a matlab expression
 % 'this' refers to the initial object/structure
  self  = this; % in case we use a Pythonic syntax
  value = value(8:end);
  try
    value = eval(value);
  catch ME
    try % try again without return value
      value = evalc(value);
    catch ME
      disp([ mfilename ': WARNING: evaluating matlab: ' value ])
      value = getReport(ME);
    end
  end
  

