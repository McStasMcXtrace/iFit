function v = fix(a)
% FIX    Round towards zero.
%   FIX(X) rounds the elements of X to the nearest integers towards zero.
%
% Example: s=iData(1:10)+0.1; all(double(fix(s)) == 1:10)
%
% Version: $Date$ $Version$ $Author$
% See also iData, iData/floor, iData/ceil, iData/round

v = unary(a, 'fix');
