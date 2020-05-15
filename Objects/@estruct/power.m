function c = power(a,b)
% .^   Object power.
%   A.^B is A to the B power.
%
%   C = POWER(A,B) is called for the syntax 'A .^ B'.
%
% Example: a=estruct(peaks); b=a^2; max(b) == max(a)^2
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/times, estruct/rdivide, estruct/mpower

if nargin == 1,
  b = a;
end
c = binary(a, b, 'power');
