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
%       The rank 0 corresponds with Signal/Monitor
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
  % travel through indexed references
  v = subsref_single(v, S(index), a);
end % for

