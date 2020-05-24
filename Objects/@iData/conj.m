function a = conj(a)
%  CONJ   Complex conjugate.
%    CONJ(X) is the complex conjugate of Signal of object X.
%    For a complex X, CONJ(X) = REAL(X) - i*IMAG(X).
%
% Example: s=iData([-i 0 i]); all(conj(s) == [i 0 -i])
% Version: $Date$ $Version$ $Author$
% See also iData, iData/transpose, iData/ctranspose, iData?imag, iData/real

a = unary(a, 'conj');


