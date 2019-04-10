function a = exp(a)
% b = exp(s) : exponential value of estruct object
%
%   @estruct/exp function to return the exponential value of data sets
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=exp(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/exp, estruct/log, estruct/log10, estruct/sqrt

a = unary(a, 'exp');

