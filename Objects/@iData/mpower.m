function c = mpower(a,b)
% ^   Object power.
%   A^B is A to the B power.
%
%   C = MPOWER(A,B) is called for the syntax 'A ^ B'.
%
% Example: a=iData(peaks); b=a^2; max(b) == max(a)^2
% Version: $Date$ $Version$ $Author$
% See also iData, iData/times, iData/rdivide, iData/power

if nargin == 1,
  b = a;
end
c = binary(a, b, 'power');
