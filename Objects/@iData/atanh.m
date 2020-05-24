function a = atanh(a)
%  ATANH  Inverse hyperbolic tangent.
%    ATANH(X) is the inverse hyperbolic tangent of the Signal of object X.
%
% Example: s=iData([-1 0 1]); all(atanh(s) == atanh([-1 0 1]))
% Version: $Date$ $Version$ $Author$
% See also iData, iData/cos, iData/acos, iData/sin, iData/asin, iData/tan, iData/atan

a = unary(a, 'atanh');

