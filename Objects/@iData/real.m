function a = real(a)
%  REAL   Complex real part.
%    REAL(X) is the real part of X.
%    See I or J to enter complex numbers.
%
% Example: s=iData(rand(5)+i*rand(5)); any(real(s),0)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/isreal, i, j, iData/imag

a = unary(a, 'real');