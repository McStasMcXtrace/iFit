function a = ceil(a)
% b = ceil(s) : upper integer round of estruct object
%
%   @estruct/ceil function to round the elements of 's' to the nearest integers
%   towards plus infinity.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=ceil(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/floor, estruct/ceil, estruct/round

a = unary(a, 'ceil');

