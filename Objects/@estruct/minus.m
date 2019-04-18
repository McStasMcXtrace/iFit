function c = minus(a,b)
% -   Minus.
%   A - B subtracts object B from A.
%   C = MINUS(A,B) is called for the syntax 'A - B' to compute the difference.
%
% Example: s=estruct(-10:10); c=s-1; c{0}(1) == -11
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/minus, estruct/plus, estruct/times, estruct/rdivide

if nargin ==1
  b=[];
end
c = binary(a, b, 'minus');