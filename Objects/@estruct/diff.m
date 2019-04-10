function a = diff(a)
% b = diff(s) : computes the difference long 1st axis of estruct object
%
%   @estruct/diff function to compute the difference along rows, that is the 
%     gradient for the 1st axis (rows).
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=diff(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/gradient, estruct/sum, estruct/trapz, estruct/jacobian

a = gradient(a, 1);

