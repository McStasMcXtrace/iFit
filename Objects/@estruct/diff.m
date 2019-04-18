function a = diff(a)
% DIFF Difference and approximate derivative along 1st axis.
%   B = DIFF(A) compute the difference along rows, that is the
%   gradient for the 1st axis (rows).
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/gradient, estruct/sum, estruct/trapz, estruct/jacobian

a = gradient(a, 1);