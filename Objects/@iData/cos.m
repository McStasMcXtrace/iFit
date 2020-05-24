function a = cos(a)
%  COS    Cosine of argument in radians.
%    COS(X) is the cosine of the elements of X. 
%
% Example: s=iData([0 pi/4 pi/3 pi/2]); all(abs(cos(s) - [1 0.707 1/2 0]) < 0.01)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/cos, iData/acos, iData/sin, iData/asin, iData/tan, iData/atan

a = unary(a, 'cos');

