function v = isnumeric(a)
% ISNUMERIC True for numeric arrays.
%    ISNUMERIC(A) returns true if A is a numeric array and false otherwise.
%
%    For example, integer and float (single and double) objects are numeric,
%    while logicals, strings, cell arrays, and structure arrays are not.
%
% Example: s=iData(-10:10); isnumeric(s)
% Version: $Date$ $Version$ $Author$
% See also iData/double, iData/isinteger, iData/isfloat

v = unary(a, 'isnumeric');
if iscell(v), v = cell2mat(v); end
v = logical(v);