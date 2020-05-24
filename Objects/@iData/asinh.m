function a = asinh(a)
%  ASINH  Inverse hyperbolic sine.
%    ASINH(X) is the inverse hyperbolic sine of the Signal of object X.
%
% Example: s=iData([-1 0 1]); all(asinh(s) == asinh([-1 0 1]))
% Version: $Date$ $Version$ $Author$
% See also iData, iData/cos, iData/acos, iData/sin, iData/asin, iData/tan, iData/atan

a = unary(a, 'asinh');

