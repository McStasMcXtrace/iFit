function c = times(a,b)
%  .*  Array multiply (times).
%    X.*Y denotes element-by-element multiplication.
%
%    C = TIMES(A,B) is called for the syntax 'A .* B'.
%
% Example: a=estruct(-10:10); c=a.*2; max(c{0}) == 20
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/minus, estruct/plus, estruct/times, estruct/rdivide, estruct/power
if nargin ==1
  b=[];
end

c = binary(a, b, 'times');

