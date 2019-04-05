function s = setalias(self, varargin)
% SETALIAS set object aliases (char properties)
%   SETALIAS(s, 'PropertyName','PropertyValue') sets the definition of the specified
%   property. The PropertyValue must be a char in order to be interpreted as an alias.
%   SETALIAS(s, ...) is equivalent to SET(s, ..., 'alias')
%
%   SETALIAS is similar to SET except when the alias value is given as a string/char.
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
%   SETALIAS(S) displays all object properties defines as strings.
%
%   SETALIAS(s, 'P1', 'V1', 'P2','V2' ...) sets multiple aliases.
%
%   SETALIAS(s, 'PropertyName')
%   SETALIAS(s, 'PropertyName','') remove PropertyName from object.
%
% Example: s=estruct; set(s, 'a', 1); setalias(s, 'b.c','a'); get(s, 'b.c') == 1 && strcmp(getalias(s, 'b.c'),'a')
% Version: $Date$ $
% See also estruct, fieldnames, findfield, isfield, set, get, getalias, setalias, 
%   getaxis, setaxis

if nargin == 1
  s = findfield(self);
  return
elseif nargin ==2
  varargin{end+1}='';
end

s = set(self, varargin{:}, 'alias');

