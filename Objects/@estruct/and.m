function c = and(a,b)
% &  Logical AND (and).
%    A & B performs a logical AND of objects A and B and returns an array
%    containing elements set to either logical 1 (TRUE) or logical 0
%    (FALSE).  An element of the output object Signal is set to 1 if both
%    input objects contain a non-zero element at that same Signal location.
%    Otherwise, that element is set to 0.
%
%    C = AND(A,B) is called for the syntax 'A & B'.
%
% Example: a=estruct(1:10); all(a & 1, 0)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/or, estruct/xor
if nargin ==1
  b=[];
end

c = binary(a, b, 'and');

