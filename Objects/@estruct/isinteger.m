function v = isinterger(a)
%  ISINTEGER True for objects of integer data type.
%    ISINTEGER(A) returns true if A is an object Signal of integer data
%    type and false otherwise.
%
%    The 8 integer data types in MATLAB are int8, uint8, int16, uint16,
%    int32, uint32, int64 and uint64.
%
% Example: s=estruct(uint8(-10:10)); isinteger(s)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar,
%          estruct/isvector, estruct/issparse

v = unary(a, 'isinteger');
if iscell(v), v = cell2mat(v); end
v = logical(v);