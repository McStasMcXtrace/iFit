function v = isinf(a)
%  ISINF  True for infinite elements.
%    ISINF(X) returns an object with Signal as 1's where the
%    Signal of X is +Inf or -Inf and 0's where it is not.
%
%    To remove nan's and inf's values use: fill(s)
%
% Example: s=estruct([pi NaN Inf -Inf]); all(isinf(s) == [0 0 1 1])
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/isfinite, estruct/isnan, estruct/fill
v = unary(a, 'isinf');