function a = sin(a)
% SIN    Sine of argument in radians.
%   SIN(X) is the sine of the elements of X.
%
% Example: s=estruct([0 pi/4 pi/6 pi/2]); all(abs(sin(s) - [0 0.7071 1/2 1]) < 0.01)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'sin');

