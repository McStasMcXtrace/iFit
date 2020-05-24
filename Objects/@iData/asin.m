function b = asin(a)
%  ASIN   Inverse sine, result in radians.
%    ASIN(X) is the arcsine of the Signal of X. Complex
%    results are obtained if ABS(x) > 1.0 for some element.
%
% Example: s=iData([-1 0 1]); all(asin(s)==[-1/2 0 1/2]*pi)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/cos, iData/acos, iData/sin, iData/asin, iData/tan, iData/atan

b=unary(a, 'asin');

