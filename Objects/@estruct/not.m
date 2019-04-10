function a = not(a)
% b = not(s) : computes the logical 'not' of estruct object
%
%   @estruct/not function to compute the logical 'not' of data sets, 
%     i.e. returns 0 where estruct Signal is not zero, and 1 otherwise.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=not(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/or, estruct/and

a = unary(a, 'not');

