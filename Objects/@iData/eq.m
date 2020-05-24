function c = eq(a,b)
%  ==  Equal (eq).
%    A == B does element by element comparisons between A and B
%    and returns an object of the same size with Signal set to logical 1
%    where the relation is true and Signal set to logical 0 where it is
%    not.
%    When comparing two iData objects, the monitor weighting is applied.
%
%    C = EQ(A,B) is called for the syntax 'A == B'.
%
% Example: s=iData(-10:10); any(s == 1)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/find, iData/gt, iData/lt, iData/ge, iData/le, iData/ne, iData/eq

if nargin ==1
  b=[];
end
c = binary(a, b, 'eq');

