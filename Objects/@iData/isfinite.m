function v = isfinite(a)
%  ISFINITE True for finite elements in Signal.
%    ISFINITE(X) returns an object with Signal as 1's where
%    the Signal of X is finite (not nan's or inf's) and 0's where
%    it is not.
%
%    For any X, exactly one of ISFINITE(X), ISINF(X), or ISNAN(X)
%    is 1 for each element of the object Signal.
%
%    To remove nan's and inf's values use: fill(s)
%
% Example: s=iData([pi NaN Inf -Inf]); all(isfinite(s) == [1 0 0 0])
% Version: $Date$ $Version$ $Author$
% See also iData, iData/isinf, iData/isnan, iData/fill

v = unary(a, 'isfinite');