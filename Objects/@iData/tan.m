function a = tan(a)
% TAN    Tangent of argument in radians.
%   TAN(X) is the tangent of the elements of X, i.e. sin(X)/cos(X).
%
% Example: s=iData([0 pi/3 pi/4 pi/6]); all(abs(tan(s) - [0 1.732 1 0.577]) < 0.01)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/cos, iData/acos, iData/sin, iData/asin, iData/tan, iData/atan

a = unary(a, 'tan');

