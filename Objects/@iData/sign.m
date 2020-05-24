function a = sign(a)
% SIGN   Signum function.
%   For each Signal element of X, SIGN(X) returns 1 if the element
%   is greater than zero, 0 if it equals zero and -1 if it is
%   less than zero.  For the nonzero elements of complex X,
%   SIGN(X) = X ./ ABS(X).
%
% Example: a=iData(peaks); b=sign(a); all(b(:)==1 | b(:)==-1)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/abs

a = unary(a, 'sign');

