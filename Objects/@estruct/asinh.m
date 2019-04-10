function a = asinh(a)
% b = asinh(s) : computes the inverse hyperbolic sine of estruct object
%
%   @estruct/asinh function to compute the inverse hyperbolic sine of data sets.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=asinh(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'asinh');

