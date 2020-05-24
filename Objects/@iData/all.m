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
%    ALL(X, 0) does the same as above, but returns a single scalar per object.
%
% Example: s=iData(-10:10); any(s<0) && all(abs(s)>=0)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/any

if nargin > 1 && isequal(varargin{1},0)
  flag_scalar = true;
  varargin(1)=[];
else
  flag_scalar = false;
end

b = unary(a, 'all', varargin{:});

if flag_scalar
  if iscell(b)
    for index=1:numel(b)
      this = b{index};
      b{index} = all(this(:));
    end
  else
    b = all(b(:));
  end
end

