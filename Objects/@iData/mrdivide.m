function c = mrdivide(a,b)
% /   Slash or right matrix divide.
%    A/B is the object division of B into A. 
%    The syntax A/B is the same as A./B
%
%    C = MRDIVIDE(A,B) is called for the syntax 'A / B'.
%
% Example: a=iData(peaks); b=a/a; all(b(:)==1)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/minus, iData/plus, iData/times, iData/rdivide, iData/power

c = rdivide(a,b);


