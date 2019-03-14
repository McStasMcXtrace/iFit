function v = subsref(a,S)
% a = subsref(a,index) indexed reference
%
%   @estruct/subsref: get indexed property 
%   such as:
%     a(1:2,3)
%       Get indexed elements in array
%
%     a.field
%     a.('field')
%       Get field in structure, same as a.b. Follow links.
%       This syntax gets final values (travel through aliases/links).
%       Equivalent to get(s, 'field')
%
%     a{axis_rank}
%       Get an axis of given rank.
%       The rank 0 corresponds with Sinal/Monitor
%       Equivalent to getaxis(s, axis_rank)
%
%     a('field')
%       Get complex field in structure and does NOT follow links.
%       This syntax gets 'aliases' (does not travel through links).
%       Equivalent to get(s, 'field', 'alias') and getalias(s, 'field')

%   In the alias case above, a string/char value allows to link to internal or 
%   external links, as well as evaluated expression, with the following syntax cases:
%     'field'                           a simple link to an other property 'field'
%     'field1.field2...'                a nested link to an other property
%     'file://some_file_path'           a local file URL
%     'http://some_distant_resource'    an HTTP URL (proxy settings may have to be set)
%     'https://some_distant_resource'   an HTTPS URL (proxy settings may have to be set)
%     'ftp://some_distant_resource'     an FTP URL (proxy settings may have to be set)
%     'matlab: some_expression'         some code to evaluate. 'this' refers to the object itself
%
%   File and URL can refer to compressed resources (zip, gz, tar, Z) which are 
%   extracted on-the-fly. In case the URL/file resource contains 'sections', a 
%   search token can be specified with syntax such as 'file://filename#token'.
%
% Version: $Date$
v = [];
if numel(a) > 1
  if strcmp(S(1).type,'()') && iscellstr(S(1).subs) && ~strcmp(S(1).subs{1},':')
    v = {}; % special case array('field') returns the field for each array element
    for index=1:numel(a)
      try
        v{index} = subsref(a(index),S);
      catch
        v{index} = [];
      end
    end
    v = reshape(v, size(a));
  else
    v = builtin('subsref', a, S);
  end
  return
end

v = a;
for index=1:numel(S)
  default = true;
  switch S(index).type
  case {'()','.'} % syntax: a('fields') does not follow links (getalias).
                  % syntax: a.('field') follows links (get), can be a compound field.
    if ischar(S(index).subs) S(index).subs = cellstr(S(index).subs); end
    if iscellstr(S(index).subs)
      % follow links for '.' subsref, not for '()'
      v = get_single(v, S(index).subs{1}, S(index).type(1)=='.', a);
      default = false;
      if S(index).type(1)=='.' && ischar(v)
        v = get_single_alias(a, v);
      end
    end
  case '{}' % syntax: a{axis_rank} get axis value/alias (getaxis)
    if isa(v, 'estruct') && numel(S(index).subs{1}) == 1 % scalar numeric or char
      v = getaxis(v, S(index).subs{1});        
      default = false;
    end
  end
  if default  % other cases
    v = builtin('subsref', v, S(index));
  end
end % for

% ------------------------------------------------------------------------------
function v = get_single(s, field, follow, s0)
  % get_single get a single field, recursively, and can follow char links
  % when follow is true, the existing field value is checked for further link
  %   then the initial structure s0 is set again.

  if nargin < 3, follow=true; end
  
  % use builtin subsref for the whole path when 'not follow'
  if ~follow && any(field == '.')
    subs=textscan(field,'%s','Delimiter','.'); subs=subs{1};
    typs=cell(size(subs)); [typs{:}] = deal('.');
    S = struct('type',typs, 'subs', subs);
    v = builtin('subsref', s, S);
    return
  else
    % when 'follow', we need to access values recursively and check for possible 'alias' (char)
    % cut the field into pieces with '.' as separator
    if any(field == '.')
      field = textscan(field,'%s','Delimiter','.'); field=field{1};
    else
      field = cellstr(field);
    end
    % now handle each level and test for char/alias value
    v = s;
    for index=1:numel(field)
      % subsref is faster than getfield which itself calls subsref
      v = builtin('subsref',v, struct('type','.','subs', field{index})); % get true value/alias (no follow)
      if follow && ischar(v)
        v = get_single_alias(s0, v); % access a link/alias in initial object/structure
      end
    end
  end
  
% ----------------------------------------------------------------------------
function v = get_single_alias(s, v)
  if isfield(s, v) % a link to a valid field
    v = builtin('subsref',s, struct('type','.','subs', v)); % get true value/alias (no follow)
  elseif strncmp(v, 'matlab',6) && numel(v) > 8 % URL matlab:
    v = get_single_eval(s, v);
  elseif exist(v, 'file') || any(strcmp(tok, {'http' 'https' 'ftp' 'file'})) % URL http https ftp file: 
    try
      v = iLoad(v);
    end
  end
  
% ----------------------------------------------------------------------------
function value = get_single_eval(this, value)
 % get_single_eval a sandbox to valuate a matlab expression
 % 'this' refers to the initial object/structure
  try
    if iscellstr(value), value = sprintf('%s\n', value{:}); end
    value = eval(value(8:end));
  catch ME
    disp([ mfilename ': WARNING: evaluating: ' value ])
    disp(getReport(ME))
  end
  

