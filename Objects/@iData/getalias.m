function s = getalias(self, varargin)
% GETALIAS Get object aliases (definition).
%   GETALIAS(s, 'PropertyName') returns the definition of the specified property.
%   GETALIAS(s, ...) is equivalent to GET(s, ..., 'alias')
%
%   GETALIAS is similar to GET except when the alias value is given as a string/char.
%   In this case, the value allows to link to internal or external links,
%   as well as evaluated expression, with the following syntax cases:
%     'field'                           a simple internal link to an other property 'field'
%     'field1.field2...'                a nested internal link to an other property
%     'file://some_file_path'           a local file URL
%     'http://some_distant_resource'    an HTTP URL (proxy settings may need to be set)
%     'https://some_distant_resource'   an HTTPS URL (proxy settings may need to be set)
%     'ftp://some_distant_resource'     an FTP URL (proxy settings may need to be set)
%     'matlab: some_expression'         some code to evaluate. 'this' refers to the object itself
%
%   File and URL can refer to compressed resources (zip, gz, tar, Z) which are
%   extracted on-the-fly. In case the URL/file resource contains 'sections' a search token
%   can be specified with syntax such as 'file://filename#token'.
%
%   GETALIAS(S) displays all object properties defined as strings.
%
%   GETALIAS(s, 'PropertyName1', 'PropertyName2', ...) return multiple aliases.
%
% Example: s=iData; set(s, 'a', 1); setalias(s, 'b.c','a'); get(s, 'b.c') == 1 && strcmp(getalias(s, 'b.c'),'a')
% Version: $Date$ $Version$ $Author$
% See also iData, fieldnames, findfield, isfield, set, get, getalias, setalias,
%   getaxis, setaxis

if nargin == 1
  [~,s] = findstr(self);
  return
end

s = get(self, varargin{:}, 'alias');
