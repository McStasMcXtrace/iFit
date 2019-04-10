function a = sinh(a)
% b = sinh(s) : computes the hyperbolic sine of estruct object
%
%   @estruct/sinh function to compute the hyperbolic sine of data sets.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=sinh(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'sinh');

