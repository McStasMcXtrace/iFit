function v = round(a)
% ROUND  Round towards nearest integer.
%   ROUND(X) rounds the elements of X to the nearest integers.
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/floor, estruct/ceil, estruct/round

v = unary(a, 'round');

