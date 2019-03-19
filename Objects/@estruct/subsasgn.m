function a = subsasgn(a,S,val)
% a = subsasgn(a,index,b) indexed assignement
%
%   @estruct/subsasgn: indexed property assignement 
%   such as:
%     a(1:2,3) =val
%       Set indexed elements in array
%
%     a.field     =val
%     a.('field') =val
%       Set field in structure, same as a.b. Follow links.
%       This syntax sets final values (travel through aliases/links).
%       Equivalent to set(s, 'field', val)
%
%     a{axis_rank}=val
%       Set an axis of given rank to a value or alias/link (when given as string/char).
%       The rank 0 corresponds with Signal/Monitor
%       Equivalent to setaxis(s, axis_rank, value)
%
%     a('field')  =val
%       Set complex field in structure and does NOT follow links.
%       This syntax sets 'aliases' (does not travel through links) when the 
%         value is a string/char.
%       Equivalent to set(s, 'field', val, 'alias') and setalias(s, 'field', val)

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

% we use a recursive approach until we reach the last level for assigment
% then propagate back the update to the upper levels
a = subsasgn_recursive(a, S, val);

% reset cache (as we have changed the object: fields, values, ...)
a.Private.cache.findfield = [];
a.ModificationDate = clock;

