function a = uminus(a)
% -  Unary minus.
%   -A negates object A.
%
%   B = UMINUS(A) is called for the syntax '-A'.
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/uminus, estruct/abs, estruct/real, estruct/imag, estruct/uplus

a = unary(a, 'uminus');