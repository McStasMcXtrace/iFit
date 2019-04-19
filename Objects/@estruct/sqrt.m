function a = sqrt(a)
%  SQRT   Square root.
%   SQRT(X) is the square root of the elements of X. Complex 
%   results are produced if X is not positive.
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sqrt, estruct/power

a = unary(a, 'sqrt');

