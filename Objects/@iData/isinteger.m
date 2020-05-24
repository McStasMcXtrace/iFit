function v = isinterger(a)
%  ISINTEGER True for objects of integer data type.
%    ISINTEGER(A) returns true if A is an object Signal of integer data
%    type and false otherwise.
%
%    The 8 integer data types in MATLAB are int8, uint8, int16, uint16,
%    int32, uint32, int64 and uint64.
%
% Example: s=iData(uint8(-10:10)); isinteger(s)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/isreal, iData/isfinite, iData/isnan,
%          iData/isinf, iData/isfloat, iData/isinterger,
%          iData/isnumeric, iData/islogical, iData/isscalar,
%          iData/isvector, iData/issparse

v = unary(a, 'isinteger');
if iscell(v), v = cell2mat(v); end
v = logical(v);