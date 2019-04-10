function a = imag(a)
% b = imag(s) : imaginary part of estruct object
%
%   @estruct/imag function to return the imaginary part of data sets.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=imag(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/uminus, estruct/abs, estruct/real, estruct/imag, estruct/uplus

a = unary(a, 'imag');

