function a = floor(a)
% b = floor(s) : lower integer round of estruct object
%
%   @estruct/floor function to round the elements of 's' to the nearest integers
%   towards minus infinity.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=floor(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/floor, estruct/ceil, estruct/round

a = unary(a, 'floor');

