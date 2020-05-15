function c = rdivide(a,b)
% ./  Right object divide.
%    A./B denotes element-by-element division (object ratio). 
%    The syntax A/B is the same as A./B
%
%    C = RDIVIDE(A,B) is called for the syntax 'A ./ B
%
% Example: a=estruct(peaks); b=a./a; all(b(:)==1)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/minus, estruct/plus, estruct/times, estruct/rdivide

c = binary(a, b, 'rdivide');

