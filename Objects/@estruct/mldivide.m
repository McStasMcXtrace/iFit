function c = mldivide(a,b)
% \  Combine (left object divide, mldivide)
%   C = MLDIVIDE(A,B) is called for the syntax 'A \ B'. It is equivalent
%   to COMBINE, so that A\B = MLDIVIDE(A,B) = COMBINE(A,B)
%
% Version: $Date$ $Version$ $Author$
% See also iData, iData/combine
if nargin ==1
  b=[];
end
c = combine(a,b);