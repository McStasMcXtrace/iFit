function b = var(a, varargin)
%  VAR Variance of object Signal.
%    S = VAR(X) returns the variance of the Signal of the object X. If
%    X is a N-D object, var(X) operates along the first
%    non-singleton dimension.
%
%    VAR normalizes Y by N-1 if N>1, where N is the sample size.  This is
%    an unbiased estimator of the variance of the population from which X is
%    drawn, as long as X consists of independent, identically distributed
%    samples. For N=1, Y is normalized by N. 
%
%    Y = VAR(X,1) normalizes by N and produces the second moment of the
%    sample about its mean.
% 
%    Y = VAR(X,W) computes the variance using the weight vector W.  The
%    length of W must equal the length of the dimension over which VAR
%    operates, and its elements must be nonnegative.  VAR normalizes W to
%    sum to one.
% 
%    Y = VAR(X,W,DIM) takes the variance along the dimension DIM of X.  Pass
%    in 0 for W to use the default normalization by N-1, or 1 to use N.
%
%    S = var(X, 0) does the same as above, but returns the total variance per object.
% 
%    The variance is the square of the standard deviation (STD).
%
% Example: s=iData(-10:10); var(s,0) == 38.5
% Version: $Date$ $Version$ $Author$
% See also iData, iData/uminus, iData/abs, iData/real, iData/imag, iData/uplus

b = unary(a, mfilename, varargin{:});

if nargin > 1 && isequal(varargin{1},0)
  if iscell(b)
    b = cell2mat(b);
  end
end

