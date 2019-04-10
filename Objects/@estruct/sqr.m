function c = sqr(a)
% c = sqr(a) : computes the square of estruct objects
%
%   @estruct/sqr function to compute the square of data sets
%
% input:  a: object or array (estruct)
% output: c: object or array (estruct)
% ex:     c=sqr(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/times, estruct/power

c = binary(a, 2, 'power');

