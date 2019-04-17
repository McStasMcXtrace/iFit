function c = xor(a,b)
%  XOR Logical EXCLUSIVE OR.
%    XOR(S,T) is the logical symmetric difference of objects S and T.
%    The result is logical 1 (TRUE) where either S or T, but not both, is
%    nonzero.  The result is logical 0 (FALSE) where S and T are both zero
%    or nonzero.
%
% Example: s=estruct(-10:10); all(xor(s, ~s), 0)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/and, estruct/or
if nargin ==1
  b=[];
end

c = binary(a, b, 'xor');

