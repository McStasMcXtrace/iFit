function a = acos(a)
% b = cos(s) : computes the cosine of estruct object
%
%   @estruct/cos function to compute the cosine of data sets (using radians).
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=cos(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'cos');

