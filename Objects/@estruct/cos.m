function a = cos(a)
%  COS    Cosine of argument in radians.
%    COS(X) is the cosine of the elements of X. 
%
% Example: s=estruct([0 pi/4 pi/3 pi/2]); all(abs(cos(s) - [1 0.707 1/2 0]) < 0.01)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'cos');

