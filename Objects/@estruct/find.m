function b = find(a, varargin)
%  FIND   Find indices of nonzero Signal elements.
%    I = FIND(X) returns the linear indices corresponding to 
%    the nonzero entries of the array X.  X may be a logical expression. 
%    Use IND2SUB(SIZE(X),I) to calculate multiple subscripts from 
%    the linear indices I.
%    Use X.Signal(I) to get the corresponding Signal elements.
%  
%    I = FIND(X,K) returns at most the first K indices corresponding to 
%    the nonzero entries of the object Signal X.  K must be a positive integer, 
%    but can be of any numeric type.
% 
%    I = FIND(X,K,'first') is the same as I = FIND(X,K).
% 
%    I = FIND(X,K,'last') returns at most the last K indices corresponding 
%    to the nonzero entries of the array X.
%
% Example: s=estruct(-10:10); numel(find(s>0)) == 10
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/uminus, estruct/abs, estruct/real, estruct/imag, estruct/uplus

b = unary(a, 'find', varargin{:});

