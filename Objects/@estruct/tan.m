function a = tan(a)
% b = tan(s) : computes the tangent of estruct object
%
%   @estruct/atan function to compute the tangent of data sets (using radians).
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=tan(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'tan');

