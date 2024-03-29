function v = subsref_single(v, S, a)
% subsref_single single level subsref
%   SUBSREF_SINGLE(v, S, a) gets v.(S) from root object a. Default with v=a.
  if nargin<3, a=v; end
  if ischar(S), S=struct('type','.','subs', S); end
  if ~isstruct(S),  error([ mfilename ': invalid reference (2nd arg) expect struct, is ' class(S) ]); end
  if numel(S) ~= 1, error([ mfilename ': only works with a single level reference' ]); end

  if isa(a, 'iData') && a.verbose > 2
    if ischar(S.subs) || iscellstr(S.subs)
      disp([ mfilename ': DEBUG: get object "' S.type char(S.subs) '"' ])
    else
      disp([ mfilename ': DEBUG: get object "' S.type '"' ])
      disp(S.subs)
    end
  end

  default = true;
  switch S.type
  case {'()','.'} % syntax: a('fields') does not follow links (getalias).
                  % syntax: a.('field') follows links (get), can be a compound field.
    if ischar(S.subs) S.subs = cellstr(S.subs); end
    if iscellstr(S.subs) && isscalar(S.subs) && (~isobject(v) || ~ismethod(v, S.subs{1}))
      % follow links for '.' subsref, not for '()'
      v = get_single(v, S.subs{1}, S.type(1)=='.', a);  % which handles aliases in 'a'
      default = false;
    elseif iscell(S.subs)
      if any(cellfun(@(c)isa(c,'iData'), S.subs))
        for index=1:numel(S.subs)
          if isa(S.subs{index},'iData')
            S.subs{index} = subsindex(S.subs{index});
          end
        end
      end
      if isa(v,'iData')
        % v = subsref_single(v,'Signal');
        v = subsref_iData(v, S);
        return
      end
    end
  case '{}' % syntax: a{axis_rank} get axis value/alias (getaxis)
    if isa(v, 'iData') && numel(S.subs{1}) == 1 % scalar numeric or char
      try
        v = getaxis(v, S.subs{1});
        default = false;
      end
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
  if iscellstr(field), field = char(field); end

  if isfield(s0, 'verbose') && s0.verbose > 2
    disp([ mfilename ': DEBUG: get_single ' field ' ' num2str(follow) ]);
  end

  % cut the field into pieces with '.' as separator
  if any(field == '.')
    field = textscan(field,'%s','Delimiter','.'); field=field{1};
    typs=cell(size(field)); [typs{:}] = deal('.');
    S = struct('type',typs, 'subs', field);
  elseif strcmp(field, ':') && ~follow % special case for (:)
    v = reshape(s, [ prod(size(s)) 1 ]);
    return
  else % access a simple field name
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
      failed = true;
      try
        v = builtin('subsref',v, S(index)); % get true value/alias (no follow)
      catch
        f = findfield(s0, S(index).subs, 'exact case');
        if ischar(f) || numel(f) == 1
          v = get_single(v, f, S.type(1)=='.', s0);
        else
          error([ mfilename ': No appropriate method, property, or unique field "' char(field) '" for class iData.' ]);
        end
      end
      % when 'follow', we need to access values iteratively and check for possible 'alias' (char)
      if follow && ischar(v) && ~strcmp(S(index).subs,'Source')
        if any(v == '.') && all(isletter(v) | v == '_' | v == '.' | isstrprop(v, 'digit')) % a compound link
          try % this may fail when e.g. field is a file name with extension -> has a dot.
            v = get_single(s, v, follow, s0);
          end
        else
          v = get_single_alias(s0, v); % access a link/alias in initial object/structure
        end
      end
    end
  end

% ----------------------------------------------------------------------------
function v = get_single_alias(s, v)
  if ~ischar(v), return; end
  if size(v,1) > 1, v=v(:)'; end % make it a row
  if isfield(s, v) % a link to a valid field in root object (single level)
    v = builtin('subsref',s, struct('type','.','subs', v)); % get true value/alias (no follow)
  elseif strcmp(v, 'matlab: sqrt(this.Signal)') % for default Error (no eval)
    v = sqrt(abs(double(get(s,'Signal'))));
  elseif strncmp(v, 'matlab',6) && numel(v) > 8 % URL matlab:
    v = get_single_eval(s, v);
  elseif any(strcmp(strtok(v, ':'), {'http' 'https' 'ftp' 'file'})) % URL http https ftp file:
    try
      v = iLoad(v);
    end
  end

% ----------------------------------------------------------------------------
function value = get_single_eval(this, value)
 % GET_SINGLE_EVAL a sandbox to evaluate a matlab expression
 % 'this' refers to the initial object/structure
  self  = this; % in case we use a Pythonic syntax
  Data  = self.Data;
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

% ----------------------------------------------------------------------------
function s = subsref_iData(a, S)

  s = copyobj(a);

  % perform rebinning on Signal
  set(s, 'Signal', subsref(get(a, 'Signal'), S));

  % perform rebinning on Error
  try
    set(s, 'Error', subsref(get(a, 'Error'), S));
  end

  % perform rebinning on Monitor
  try
    set(s, 'Monitor', subsref(get(a, 'Monitor'), S));
  end

  % perform rebinning on axes
  myisvector = @(c)length(c) == numel(c);
  for index=1:length(S.subs)
    x = getaxis(a, index);
    try
      if myisvector(x),
        if ~strcmp(S.subs{index},':'), x = x(S.subs{index}); end
      else
        x = subsref(x,S);
      end
      setaxis(s, index, x);
    end
  end
  s.Private.cache.size = [];
  s.Private.cache.std_c= [];
  s.Private.cache.std_w= [];
  axescheck(s);
