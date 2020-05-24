function v = round(a)
% ROUND  Round towards nearest integer.
%   ROUND(X) rounds the elements of X to the nearest integers.
%
% Example: x=-5:5; a=iData(x+rand/3); all(round(a) == x)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/floor, iData/ceil, iData/round

v = unary(a, 'round');

