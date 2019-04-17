function v = isfloat(a)
%  ISFLOAT True for floating point object, both single and double.
%    ISFLOAT(A) returns true if A is a floating point object and false otherwise.
%
%    Single and double are the only floating point data types in MATLAB.
%
% Example: s=estruct(-10:10); isfloat(s)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sign, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar,
%          estruct/isvector, estruct/issparse

v = unary(a, 'isfloat');
if iscell(v), v=cell2mat(v); end
