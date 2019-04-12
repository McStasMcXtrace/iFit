function a = atan(a)
%  ATAN   Inverse tangent, result in radians.
%    ATAN(X) is the arctangent of the Signal of object X.
%
% Example: s=estruct([-1 0 1]); all(atan(s)==[-1 0 1]*pi/4)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'atan');

