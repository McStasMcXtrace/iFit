function a = log10(a)
% b = log10(s) : base 10 logarithm value of estruct object
%
%   @estruct/log10 function to return the base 10 logarithm value of data sets
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=log10(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/exp, estruct/log, estruct/log10, estruct/sqrt

a = unary(a, 'log10');

