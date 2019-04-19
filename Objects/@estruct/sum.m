function b = sum(a, varargin)
%  SUM Sum of in object Signal.
%    S = SUM(X) is the sum of the Signal of the object X. If
%    X is a N-D object, SUM(X) operates along the first
%    non-singleton dimension.
%    If X is floating point, that is double or single, S is
%    accumulated natively, that is in the same class as X,
%    and S has the same class as X. If X is not floating point,
%    S is accumulated in double and S has class double.
%
%    S = SUM(X,DIM) sums along the dimension DIM. 
%
%    S = SUM(X,'double') and S = SUM(X,DIM,'double') accumulate
%    S in double and S has class double, even if X is single.
%
%    S = SUM(X,'native') and S = SUM(X,DIM,'native') accumulate
%    S natively and S has the same class as X.
%
%    S = SUM(X, 0) does the same as above, but returns the total sum per object.
%
% Example: s=estruct(-10:10); sum(s,0) == 0
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/uminus, estruct/abs, estruct/real, estruct/imag, estruct/uplus

if nargin > 1 && isequal(varargin{1},0)
  flag_scalar = true;
  varargin(1)=[];
else
  flag_scalar = false;
end

b = unary(a, 'sum', varargin{:});

if flag_scalar
  if iscell(b)
    for index=1:numel(b)
      this = b{index};
      b{index} = sum(this(:));
    end
  else
    b = sum(b(:));
  end
end

