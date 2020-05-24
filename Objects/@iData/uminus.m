function a = uminus(a)
% -  Unary minus.
%   -A negates object A.
%
%   B = UMINUS(A) is called for the syntax '-A'.
%
% Example: a=iData(-5:-1); all(-a == 5:-1:1)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/uminus, iData/abs, iData/real, iData/imag, iData/uplus

a = unary(a, 'uminus');
