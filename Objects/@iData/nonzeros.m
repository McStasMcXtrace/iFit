function b = nonzeros(a)
%   NONZEROS Nonzero matrix elements.
%    NONZEROS(S) is a full column vector of the nonzero elements of S Signal.
%    This gives the values, but not the indices from e.g. find(S).
%
% Example: s=iData(-10:10); numel(nonzeros(s)) == 20
% Version: $Date$ $Version$ $Author$
% See also iData, iData/uminus, iData/abs, iData/real, iData/imag, iData/uplus

b = unary(a, 'nonzeros');

