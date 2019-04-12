function a = cosh(a)
%  COSH   Hyperbolic cosine.
%    COSH(X) is the hyperbolic cosine of the Signal of object X.
%
% Example: s=estruct(0:10); all(cosh(s) == cosh(0:10))
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'cosh');

