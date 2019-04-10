function a = log(a)
% b = log(s) : natural logarithm value of estruct object
%
%   @estruct/log function to return the natural logarithm value of data sets, ln(s)
%   Use log10 for the base 10 log.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=log(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/exp, estruct/log, estruct/log10, estruct/sqrt

a = unary(a, 'log');

