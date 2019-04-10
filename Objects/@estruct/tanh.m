function a = tanh(a)
% b = tanh(s) : computes the hyperbolic tangent of estruct object
%
%   @estruct/tanh function to compute the hyperbolic tangent of data sets.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=tanh(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'tanh');

