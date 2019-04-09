function b = any(a, varargin)
% ANY    True if any element of a object Signal is a nonzero number or is
%    logical 1 (TRUE).  ANY ignores entries that are NaN (Not a Number).
% 
%    For vector objects (1D), ANY(V) returns logical 1 (TRUE) if any of the 
%    Signal of the object is a nonzero number or is logical 1 (TRUE).
%    Otherwise it returns logical 0 (FALSE).  For matrix (2D), ANY(X) 
%    operates on the columns of X Signal, returning a row vector of logical 1's 
%    and 0's.  For multi-dimensional objects, ANY(X) operates on the 
%    first non-singleton dimension.
% 
%    ANY(X,DIM) works down the dimension DIM.  For example, ANY(X,1)
%    works down the first dimension (the rows) of X.
%
% Example: s=estruct(-10:10); any(s.Signal< 0) && all(double(abs(s))>=0)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/all

b = unary(a, 'any', varargin{:});

