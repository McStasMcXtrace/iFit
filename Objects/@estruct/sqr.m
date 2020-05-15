function c = sqr(a)
% SQR computes the square of object.
%   SQR(X) has its Signal squared, and is equivalent to POWER(X,2).
%
% Example: s=estruct(0:10); all(sqr(s) == (0:10).^2)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/times, estruct/power

c = binary(a, 2, 'power');

