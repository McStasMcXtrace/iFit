function b = acos(a)
%  ACOS   Inverse cosine, result in radians.
%    ACOS(X) is the arccosine of the elements of X. Complex
%    results are obtained if ABS(x) > 1.0 for some element.
%
% Example: s=iData([-1 0 1]); all(acos(s)==[1 1/2 0]*pi)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/cos, iData/acos, iData/sin, iData/asin, iData/tan, iData/atan

b = unary(a, 'acos');

