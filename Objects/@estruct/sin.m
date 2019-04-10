function a = sin(a)
% b = sin(s) : computes the sine of estruct object
%
%   @estruct/acos function to compute the sine of data sets (using radians).
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=sin(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'sin');

