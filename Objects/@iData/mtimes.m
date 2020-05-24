function c = mtimes(a,b)
% *  Object multiply (mtimes).
%   X*Y denotes element-by-element multiplication. 
%   The syntax X*Y is the same as X.*Y
%
%   C = MTIMES(A,B) is called for the syntax 'A * B'.
%
% Example: a=iData(-10:10); c=a.*2; max(c{0}) == 20
% Version: $Date$ $Version$ $Author$
% See also iData, iData/minus, iData/plus, iData/times, iData/rdivide, iData/power

c = times(a,b);
