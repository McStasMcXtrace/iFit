function a = uplus(a)
% b = uplus(s) : makes a copy of estruct object
%
%   @estruct/uplus function to return a duplicate of data sets.
%   b=+a creates a new estruct object with same content as 'a', but different Tag/ID and Date.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=uplus(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/uminus, estruct/abs, estruct/real, estruct/imag, estruct/uplus

a = unary(a, 'uplus');

