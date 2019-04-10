function a = del2(a)
% b = del2(s) : computes the Discrete Laplacian of estruct object
%
%   @estruct/del2 function to compute the Discrete Laplacian of data sets.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=del2(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/gradient, del2, gradient, estruct/jacobian

% make sure axes are regularly binned
%a = interp(a);

a = unary(a, 'del2');

