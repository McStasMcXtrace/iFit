function a = round(a)
% b = round(s) : upper integer round of estruct object
%
%   @estruct/round function to round the elements of 's' to the nearest integers.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=round(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/floor, estruct/ceil, estruct/round

a = unary(a, 'round');

