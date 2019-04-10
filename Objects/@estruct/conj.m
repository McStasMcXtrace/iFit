function a = conj(a)
% b = conj(s) : conjugate of estruct object
%
%   @estruct/conj function to return conjugate of data sets
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=conj(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/transpose, estruct/ctranspose, estruct?imag, estruct/real

a = unary(a, 'conj');


