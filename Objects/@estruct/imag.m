function a = imag(a)
%  IMAG   Complex imaginary part.
%    IMAG(X) is the imaginary part of X.
%    See I or J to enter complex numbers.
%
% Example: s=estruct(rand(5)+i*rand(5)); any(imag(s),0)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/isreal, i, j, estruct/real

a = unary(a, 'imag');