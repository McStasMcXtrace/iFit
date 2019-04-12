function a = floor(a)
%  FLOOR  Round towards minus infinity.
%    FLOOR(X) rounds the Signal of X to the nearest integers
%    towards minus infinity.
%
% Example: s=estruct([-1.2 0.6 1.5]); all(floor(s) == [-2 0 1])
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/floor, estruct/ceil, estruct/round

a = unary(a, 'floor');

