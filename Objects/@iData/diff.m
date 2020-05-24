function a = diff(a)
% DIFF Difference and approximate derivative along 1st axis.
%   B = DIFF(A) compute the difference along rows, that is the
%   gradient for the 1st axis (rows).
%
% Example: a=iData(peaks); b=diff(a); sum(b,0) < 1
% Version: $Date$ $Version$ $Author$
% See also iData, iData/gradient, iData/sum, iData/trapz, iData/jacobian

a = gradient(a, 1);
