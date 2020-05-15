function v = round(a)
% ROUND  Round towards nearest integer.
%   ROUND(X) rounds the elements of X to the nearest integers.
%
% Example: x=-5:5; a=estruct(x+rand/3); all(round(a) == x)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/floor, estruct/ceil, estruct/round

v = unary(a, 'round');

