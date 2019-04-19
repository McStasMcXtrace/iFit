function c = mtimes(a,b)
% *  Object multiply (mtimes).
%   X*Y denotes element-by-element multiplication. 
%   The syntax X*Y is the same as X.*Y
%
%   C = MTIMES(A,B) is called for the syntax 'A * B'.
%
% Example: a=estruct(-10:10); c=a.*2; max(c{0}) == 20
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/minus, estruct/plus, estruct/times, estruct/rdivide, estruct/power

c = times(a,b);
