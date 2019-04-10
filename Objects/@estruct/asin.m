function b = asin(a)
% ASIN   Inverse sine, result in radians.
%    ASIN(X) is the arcsine of the Signal of X. Complex
%    results are obtained if ABS(x) > 1.0 for some element.
%
% Example: s=estruct([-1 0 1]); all(double(asin(s))==[-1/2 0 1/2]*pi)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

b=unary(a, 'asin');

