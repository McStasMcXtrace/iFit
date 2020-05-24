function b = acosh(a)
%  ACOSH  Inverse hyperbolic cosine.
%    ACOSH(X) is the inverse hyperbolic cosine of the Signal of X.
%
% Example: s=iData([1 2 3]); all(acosh(s) == acosh([1 2 3]))
% Version: $Date$ $Version$ $Author$
% See also iData, iData/cos, iData/acos, iData/sin, iData/asin, iData/tan, iData/atan

b = unary(a, 'acosh');

