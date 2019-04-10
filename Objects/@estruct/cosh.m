function a = cosh(a)
% b = cosh(s) : computes the hyperbolic cosine of estruct object
%
%   @estruct/cosh function to compute the hyperbolic cosine of data sets.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=cosh(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'cosh');

