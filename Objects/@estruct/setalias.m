function s = setalias(self, varargin)
% setalias(s, AliasName, AliasValue, ...) set estruct aliases
% 
%   @estruct/setalias function to set estruct aliases.
%
%   SETALIAS is similar to SET except when the alias value is given as a string/char.
%   In this case, value allows to link to internal or external links, 
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
%   To remove an alias, use setalias(s, AliasName) or setalias(s, AliasName, '')
%   or rmfield(s, AliasName)
%
%   SETALIAS(s, ...) is equivalent to SET(s, ..., 'alias')
%
%   SETALIAS(S) displays all object properties.
%
% Example: s=estruct; set(s, 'a', 1); setalias(s, 'b.c','a'); get(s, 'b.c') == 1 && strcmp(getalias(s, 'b.c'),'a')
% 
%  Version: $Date$
%  See also estruct/getalias, estruct/get, estruct/set

if nargin == 1
  s = findfield(self);
  return
elseif nargin ==2
  varargin{end+1}='';
end

s = set(self, varargin{:}, 'alias');

