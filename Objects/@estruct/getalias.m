function s = getalias(self, varargin)
% getalias(s, AliasName, AliasValue, ...) get estruct aliases
% 
%   @estruct/setalias function to set estruct aliases.
%
%   GETALIAS is similar to GET except when the alias value is given as a string/char.
%   In this case, the value allows to link to internal or external links, 
%   as well as evaluated expression, with the following syntax cases:
%     'field'                           a simple internal link to an other property 'field'
%     'field1.field2...'                a nested internal link to an other property
%     'file://some_file_path'           a local file URL
%     'http://some_distant_resource'    an HTTP URL (proxy settings may have to be set)
%     'https://some_distant_resource'   an HTTPS URL (proxy settings may have to be set)
%     'ftp://some_distant_resource'     an FTP URL (proxy settings may have to be set)
%     'matlab: some_expression'         some code to evaluate. 'this' refers to the object itself
%
%   File and URL can refer to compressed resources (zip, gz, tar, Z) which are 
%   extracted on-the-fly. In case the URL/file resource contains 'sections' a search token
%   can be specified with syntax such as 'file://filename#token'.
%
%   GETALIAS(S) displays all object properties defines as strings.
%
% Example: s=estruct; set(s, 'a', 1); setalias(s, 'b.c','a'); get(s, 'b.c') == 1 && strcmp(getalias(s, 'b.c'),'a')
% 
%  Version: $Date$
%  See also estruct/getalias, estruct/get, estruct/set

if nargin == 1
  [~,s] = findstr(self);
  return
end

s = get(self, varargin{:}, 'alias');
