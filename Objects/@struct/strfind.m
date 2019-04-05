function [match, field] = strfind(s, varargin)
% [match, field]=strfind(s, str, option) : look for strings stored in struct
%
%   @struct/strfind function to look for strings stored in struct
%
%   [match,field] = strfind(struct, str) returns the string containg 'str' 
%     and the field name it appears in. If 'str' is set to '', the content of all
%     character fields is returned.
%   The 'option' may contain 'exact' to search for the exact occurence, and 'case'
%   to specifiy a case sensitive search.
%
% input:  s: object or array (struct)
%         str: string to search in object, or '' (char or cellstr).
%         option: 'exact' 'case' or '' (char)
% output: match: content of struct fields that contain 'str' (cellstr)
%         field: name of struct fields that contain 'str' (cellstr)
% ex:     strfind(struct,'ILL') or strfind(s,'TITLE','exact case')
%
% Version: $Date$ $Version$ $Author$
% See also struct, struct/set, struct/get, struct/findobj, struct/findfield

[match, field] = strfind(s, varargin{:});
