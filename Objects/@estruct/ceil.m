function a = ceil(a)
%  CEIL   Round towards plus infinity.
%    CEIL(X) rounds the Signal of object X to the nearest integers
%    towards infinity.
%
% Example: s=estruct([-1.2 0.6 1.5]); all(ceil(s) == [-1 1 2])
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/floor, estruct/ceil, estruct/round

a = unary(a, 'ceil');

