function b = acosh(a)
%  ACOSH  Inverse hyperbolic cosine.
%    ACOSH(X) is the inverse hyperbolic cosine of the Signal of X.
%
% Example: s=estruct([1 2 3]); all(acosh(s) == acosh([1 2 3]))
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

b = unary(a, 'acosh');

