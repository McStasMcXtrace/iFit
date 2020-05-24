function a = log10(a)
% LOG10  Common (base 10) logarithm.
%    LOG10(X) is the base 10 logarithm of the Signal of X.
%    Complex results are produced if X is not positive.
%
% Example: s=iData([10 100 1000]); all(log10(s) == [1 2 3])
% Version: $Date$ $Version$ $Author$
% See also iData, iData/exp, iData/log

a = unary(a, 'log10');