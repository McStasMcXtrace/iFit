function c = or(a,b)
%  |   Logical OR (or).
%    A | B performs a logical OR of arrays A and B and returns an array
%    containing elements set to either logical 1 (TRUE) or logical 0
%    (FALSE).  An element of the output array is set to 1 if either input
%    array contains a non-zero element at that same array location.
%    Otherwise, that element is set to 0.  A and B must have the same
%    dimensions unless one is a scalar.
%
%    C = OR(A,B) is called for the syntax 'A | B'.
%
% Example: s=estruct(-10:10); all(s | 1,0)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/and, estruct/xor
if nargin ==1
  b=[];
end

c = binary(a, b, 'or');

