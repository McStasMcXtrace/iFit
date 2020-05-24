function b = abs(a)
%  ABS Absolute value of object
%    ABS(X) return a new object with its Signal set as the absolute value of the
%    input object.
%
% Example: s=iData(-10:10); any(s< 0) && all(abs(s))>=0
% Version: $Date$ $Version$ $Author$
% See also iData, iData/uminus, iData/abs, iData/real, iData/imag, iData/uplus

b = unary(a, 'abs');