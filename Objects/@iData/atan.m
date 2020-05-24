function a = atan(a)
%  ATAN   Inverse tangent, result in radians.
%    ATAN(X) is the arctangent of the Signal of object X.
%
% Example: s=iData([-1 0 1]); all(atan(s)==[-1 0 1]*pi/4)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/cos, iData/acos, iData/sin, iData/asin, iData/tan, iData/atan

a = unary(a, 'atan');

