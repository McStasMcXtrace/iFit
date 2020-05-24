function a = tanh(a)
% TANH   Hyperbolic tangent.
%   TANH(X) is the hyperbolic tangent of X.
%
% Example: s=iData(0:10); all(tanh(s) == tanh(0:10))
% Version: $Date$ $Version$ $Author$
% See also iData, iData/cos, iData/acos, iData/sin, iData/asin, iData/tan, iData/atan

a = unary(a, 'tanh');

