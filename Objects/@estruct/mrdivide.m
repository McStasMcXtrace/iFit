function c = mrdivide(a,b)
% /   Slash or right matrix divide.
%    A/B is the object division of B into A. 
%    The syntax A/B is the same as A./B
%
%    C = MRDIVIDE(A,B) is called for the syntax 'A / B'.
%
% Example: a=estruct(peaks); b=a/a; all(b(:)==1)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/minus, estruct/plus, estruct/times, estruct/rdivide, estruct/power

c = rdivide(a,b);


