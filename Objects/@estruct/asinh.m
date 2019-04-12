function a = asinh(a)
%  ASINH  Inverse hyperbolic sine.
%    ASINH(X) is the inverse hyperbolic sine of the Signal of object X.
%
% Example: s=estruct([-1 0 1]); all(asinh(s) == asinh([-1 0 1]))
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'asinh');

