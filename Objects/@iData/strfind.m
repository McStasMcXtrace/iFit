function [match, field] = strfind(s, varargin)
% STRFIND Find a string within object.
%   STRFIND is equivalent to FINDSTR.
%
% Example: s=iData('x',1:10,'y','blah'); ischar(strfind(s, 'blah'))
% Version: $Date$ $Version$ $Author$
% See also iData, iData/findstr, set, get, findobj, findfield

[match, field] = findstr(s, varargin{:});