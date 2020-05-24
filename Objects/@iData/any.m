function b = any(a, varargin)
%  ANY    True if any element of a object Signal is a nonzero number or is
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
%    ANY(X, 0) does the same as above, but returns a single scalar per object.
%
% Example: s=iData(-10:10); any(s<0) && all(abs(s)>=0)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/all

if nargin > 1 && isequal(varargin{1},0)
  flag_scalar = true;
  varargin(1)=[];
else
  flag_scalar = false;
end

b = unary(a, 'any', varargin{:});

if flag_scalar
  if iscell(b)
    for index=1:numel(b)
      this = b{index};
      b{index} = any(this(:));
    end
  else
    b = any(b(:));
  end
end

