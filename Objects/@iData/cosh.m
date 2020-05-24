function a = cosh(a)
%  COSH   Hyperbolic cosine.
%    COSH(X) is the hyperbolic cosine of the Signal of object X.
%
% Example: s=iData(0:10); all(cosh(s) == cosh(0:10))
% Version: $Date$ $Version$ $Author$
% See also iData, iData/cos, iData/acos, iData/sin, iData/asin, iData/tan, iData/atan

a = unary(a, 'cosh');

