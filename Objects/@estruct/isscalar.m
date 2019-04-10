function v = isscalar(a)
%  ISSCALAR True if array is a scalar.
%    ISSCALAR(S) returns logical true (1) if S is a 1 x 1 matrix
%    and logical false (0) otherwise.
%
% Example: s=estruct(pi); isscalar(s) == 1
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sign, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar, 
%          estruct/isvector, estruct/issparse

v = unary(a, 'isscalar');

