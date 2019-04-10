function a = islogical(a)
% b = islogical(s) : True for logical estruct object elements
%
%   @estruct/islogical function to return true for logical elements
%   of 's', i.e. that are of type 'true (1) or false(0)'.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=islogical(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sign, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar, 
%          estruct/isvector, estruct/issparse

a = unary(a, 'islogical');
a = uint8(a);
