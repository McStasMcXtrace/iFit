function a = atanh(a)
%  ATANH  Inverse hyperbolic tangent.
%    ATANH(X) is the inverse hyperbolic tangent of the Signal of object X.
%
% Example: s=estruct([-1 0 1]); all(atanh(s) == atanh([-1 0 1]))
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'atanh');

