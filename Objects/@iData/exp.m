function a = exp(a)
%  EXP    Exponential.
%    EXP(X) is the exponential of the Signal of object X, e to the X.
%    For complex Z=X+i*Y, EXP(Z) = EXP(X)*(COS(Y)+i*SIN(Y)).
%
% Example: s=iData(0:10); all(exp(s) == exp(0:10))
% Version: $Date$ $Version$ $Author$
% See also iData, iData/exp, iData/log, iData/log10, iData/sqrt

a = unary(a, 'exp');

