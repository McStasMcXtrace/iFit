function b = isnan(a)
%  ISNAN  True for Not-a-Number (NaN).
%    ISNAN(X) returns an array that contains 1's where
%    the elements of X are NaN's and 0's where they are not.
%
%    To remove nan's and inf's values use: fill(s)
%
% Example: s=estruct([pi NaN Inf -Inf]); all(isnan(s) == [0 1 0 0])
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/isfinite, estruct/isinf, estruct/fill

b = unary(a, 'isnan');