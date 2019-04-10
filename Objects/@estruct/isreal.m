function a = isreal(a)
% b = isreal(s) : True for real estruct object elements
%
%   @estruct/isreal function to return true for real elements
%   of 's', i.e. that are not complex.
%
% input:  s: object or array (estruct)
% output: b: array (int)
% ex:     b=isreal(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sign, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar, 
%          estruct/isvector, estruct/issparse

a = unary(a, 'isreal');
a = uint8(a);
