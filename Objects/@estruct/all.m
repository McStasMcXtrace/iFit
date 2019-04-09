function b = all(a, varargin)
%  ALL    True if all elements of object Signal are nonzero.
%    For vector objects (1D), ALL(V) returns logical 1 (TRUE) if none of the Signal 
%    are zero.  Otherwise it returns logical 0 (FALSE).  For 
%    matrix objects (2D), ALL(X) operates on the columns of X, returning a row vector
%    of logical 1's and 0's. For N-D objects, ALL(X) operates on the first
%    non-singleton dimension.
% 
%    ALL(X,DIM) works down the dimension DIM.  For example, ALL(X,1)
%    works down the first dimension (the rows) of X.
%
% Example: s=estruct(-10:10); any(s.Signal< 0) && all(double(abs(s))>=0)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/any

b = unary(a, 'all', varargin{:});

