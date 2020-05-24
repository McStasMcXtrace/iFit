function [b,sigma] = sum(a, dim)
%  SUM Sum of in object Signal.
%    S = SUM(X) is the sum of the Signal of the object X. If
%    X is a N-D object, SUM(X) operates along the first
%    dimension.
%    If X is floating point, that is double or single, S is
%    accumulated natively, that is in the same class as X,
%    and S has the same class as X. If X is not floating point,
%    S is accumulated in double and S has class double.
%
%    S = SUM(X,DIM) sums along the dimension DIM. 
%
%    [S,sigma] = SUM(X, 0) does the same as above, but returns the total sum 
%    per object. In this case, a second output argument holds the error bar.
%
% Example: s=iData(-10:10); sum(s,0) == 0
% Version: $Date$ $Version$ $Author$
% See also iData, iData/uminus, iData/abs, iData/real, iData/imag, iData/uplus

if nargin < 2, dim=1; end

[b,sigma] = private_sumtrapzproj(a,dim, 'sum');

% b = unary(a, mfilename, varargin{:});

if nargin > 1 && isequal(dim,0)
  if iscell(b)
    b = cell2mat(b);
  end
end

