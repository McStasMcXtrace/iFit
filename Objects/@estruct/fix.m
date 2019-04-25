function v = fix(a)
% FIX    Round towards zero.
%   FIX(X) rounds the elements of X to the nearest integers towards zero.
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/floor, estruct/ceil, estruct/round

v = unary(a, 'fix');