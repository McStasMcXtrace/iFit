function v = isfloat(a)
%  ISFLOAT True for floating point object, both single and double.
%    ISFLOAT(A) returns true if A is a floating point object and false otherwise.
%
%    Single and double are the only floating point data types in MATLAB.
%
% Example: s=iData(-10:10); isfloat(s)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/sign, iData/isreal, iData/isfinite, iData/isnan,
%          iData/isinf, iData/isfloat, iData/isinterger,
%          iData/isnumeric, iData/islogical, iData/isscalar,
%          iData/isvector, iData/issparse

v = unary(a, 'isfloat');
if iscell(v), v=cell2mat(v); end
v = logical(v);
