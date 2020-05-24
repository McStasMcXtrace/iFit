function b = log(a)
%  LOG    Natural logarithm.
%    LOG(Z) is the natural logarithm (ln) of the elements of X.
%    Complex results are produced if Z is not positive.
%    For complex or negative Z=x+i*y, the complex logarithm is returned.
%      log(Z) = log(abs(Z)) + i*atan2(y,x)
%    Use log10 for the base 10 log.
%
% Example: s=iData([-1 -1]); all(imag(log(s)) == pi)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/exp, iData/log10

b = unary(a, 'log');