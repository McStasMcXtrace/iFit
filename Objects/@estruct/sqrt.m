function a = sqrt(a)
% b = sqrt(s) : square root value of estruct object
%
%   @estruct/sqrt function to return the square root value of data sets
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=sqrt(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sqrt, estruct/power

a = unary(a, 'sqrt');

