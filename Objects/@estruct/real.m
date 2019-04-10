function a = real(a)
% b = real(s) : real value of estruct object
%
%   @estruct/real function to return the real value of data sets.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=real(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/uminus, estruct/abs, estruct/real, estruct/imag, estruct/uplus

a = unary(a, 'real');

