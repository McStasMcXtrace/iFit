function a = atanh(a)
% b = atanh(s) : computes the inverse hyperbolic tangent of estruct object
%
%   @estruct/atanh function to compute the inverse hyperbolic tangent of data sets.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=atanh(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'atanh');

