function a = sinh(a)
% SINH   Hyperbolic sine.
%   SINH(X) is the hyperbolic sine of X.
%
% Example: s=estruct(0:10); all(sinh(s) == sinh(0:10))
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'sinh');

