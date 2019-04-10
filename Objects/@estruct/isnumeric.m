function a = isnumeric(a)
% b = isnumeric(s) : True for numeric estruct object elements
%
%   @estruct/isnumeric function to return true for numeric elements
%   of 's', i.e. that are of type double, single, integer.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=isnumeric(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sign, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar, 
%          estruct/isvector, estruct/issparse

a = unary(a, 'isnumeric');
a = uint8(a);
