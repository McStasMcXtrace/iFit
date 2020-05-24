function c = ldivide(a,b)
% .\  Combine (left object divide, ldivide)
%   C = LDIVIDE(A,B) is called for the syntax 'A .\ B'. It is equivalent
%   to COMBINE, so that A.\B = A\B = LDIVIDE(A,B) = COMBINE(A,B)
%
% Example: a=iData(peaks); a.Monitor=1; b=a.\a; max(b) == max(a)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/minus, iData/plus, iData/times, iData/mldivide
if nargin ==1
  b=[];
end
c = combine(a,b);
