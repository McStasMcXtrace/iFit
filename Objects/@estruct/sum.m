function b = sum(a, varargin)
% ABS absolute value of object
%
%   ABS(X) return a new object with its Signal set as the absolute value of the
%   input object.
%
% Example: s=estruct(-10:10); any(s< 0) && all(abs(s))>=0
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/uminus, estruct/abs, estruct/real, estruct/imag, estruct/uplus

b = unary(a, 'sum', varargin{:});

