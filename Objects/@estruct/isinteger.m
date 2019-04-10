function a = isinterger(a)
% b = isinterger(s) : True for integer estruct object elements
%
%   @estruct/isinterger function to return true for integer elements
%   of 's', i.e. that are of type 'uint' and 'int'.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=isinterger(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sign, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar, 
%          estruct/isvector, estruct/issparse

a = unary(a, 'isinteger');
a = uint8(a);
