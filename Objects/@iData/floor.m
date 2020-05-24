function a = floor(a)
%  FLOOR  Round towards minus infinity.
%    FLOOR(X) rounds the Signal of X to the nearest integers
%    towards minus infinity.
%
% Example: s=iData([-1.2 0.6 1.5]); all(floor(s) == [-2 0 1])
% Version: $Date$ $Version$ $Author$
% See also iData, iData/floor, iData/ceil, iData/round

a = unary(a, 'floor');

