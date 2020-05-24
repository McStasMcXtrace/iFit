function a = sinh(a)
% SINH   Hyperbolic sine.
%   SINH(X) is the hyperbolic sine of X.
%
% Example: s=iData(0:10); all(sinh(s) == sinh(0:10))
% Version: $Date$ $Version$ $Author$
% See also iData, iData/cos, iData/acos, iData/sin, iData/asin, iData/tan, iData/atan

a = unary(a, 'sinh');

