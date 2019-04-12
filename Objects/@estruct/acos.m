function b = acos(a)
%  ACOS   Inverse cosine, result in radians.
%    ACOS(X) is the arccosine of the elements of X. Complex
%    results are obtained if ABS(x) > 1.0 for some element.
%
% Example: s=estruct([-1 0 1]); all(acos(s)==[1 1/2 0]*pi)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

b = unary(a, 'acos');

