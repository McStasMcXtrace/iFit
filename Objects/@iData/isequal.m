function c = isequal(a,b)
%  ISEQUAL True if object Signals are numerically equal.
%    NaNs are considered equal to each other.
%    When comparing two iData objects, the monitor weighting is applied.
%
%    ISEQUAL(A,B) is 1 if the two object Signals are the same size
%    and contain the same values, and 0 otherwise.
%
% Example: a=iData([ nan -10:10 ]); isequal(a, [ nan -10:10 ])
% Version: $Date$ $Version$ $Author$
% See also iData, iData/find, iData/gt, iData/lt, iData/ge, iData/le, iData/ne, iData/eq

if nargin ==1
  b=[];
end

if exist('isequaln')
  c = binary(a, b, 'isequaln');
elseif exist('isequalwithequalnans')
  c = binary(a, b, 'isequalwithequalnans');
else
  c = binary(a, b, 'isequal');
end
if iscell(c), c = cell2mat(c); end