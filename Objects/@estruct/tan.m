function a = tan(a)
% TAN    Tangent of argument in radians.
%   TAN(X) is the tangent of the elements of X, i.e. sin(X)/cos(X).
%
% Example: s=estruct([0 pi/3 pi/4 pi/6]); all(abs(tan(s) - [0 1.732 1 0.577]) < 0.01)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'tan');

