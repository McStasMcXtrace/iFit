function a = norm(a)
% b = norm(s) : norm-2 of estruct object
%
%   @estruct/norm function to return the norm-2 of data sets
%   This function computes the norm of the object 's'.
%
% input:  s: object or array (estruct)
% output: b: norm (double)
% ex:     b=norm(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sum, estruct/trapz, norm

a = unary(a, 'norm');

