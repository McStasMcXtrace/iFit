function v = fix(a)
% FIX    Round towards zero.
%   FIX(X) rounds the elements of X to the nearest integers towards zero.
%
% Example: s=estruct(1:10)+0.1; all(double(fix(s)) == 1:10)
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/floor, estruct/ceil, estruct/round

v = unary(a, 'fix');
