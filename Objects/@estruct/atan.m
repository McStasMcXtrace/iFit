function a = atan(a)
% b = atan(s) : computes the arc tangent of estruct object
%
%   @estruct/atan function to compute the inverse tangent of data sets (in radians).
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=atan(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'atan');

