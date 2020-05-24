function c = times(a,b)
% .*  Object multiply (times).
%   X.*Y denotes element-by-element multiplication.
%   The syntax X*Y is the same as X.*Y
%
%   C = TIMES(A,B) is called for the syntax 'A .* B'.
%
% Example: a=iData(-10:10); c=a.*2; max(c{0}) == 20
% Version: $Date$ $Version$ $Author$
% See also iData, iData/minus, iData/plus, iData/times, iData/rdivide, iData/power
if nargin ==1
  b=[];
end

c = binary(a, b, 'times');
