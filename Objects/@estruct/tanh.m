function a = tanh(a)
% TANH   Hyperbolic tangent.
%   TANH(X) is the hyperbolic tangent of X.
%
% Example: s=estruct(0:10); all(tanh(s) == tanh(0:10))
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cos, estruct/acos, estruct/sin, estruct/asin, estruct/tan, estruct/atan

a = unary(a, 'tanh');

