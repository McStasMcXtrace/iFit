function a = uminus(a)
% b = uminus(s) : opposite value of estruct object
%
%   @estruct/uminus function to return the opposite value of data sets, i.e. b=-s.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=uminus(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/uminus, estruct/abs, estruct/real, estruct/imag, estruct/uplus

a = unary(a, 'uminus');

