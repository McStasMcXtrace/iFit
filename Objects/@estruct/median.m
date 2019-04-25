function b = median(a, varargin)
%  MEDIAN Median value of Signal.
%    S = MEDIAN(X) is the median value of the Signal of the object X. If
%    X is a N-D object, median(X) operates along the first
%    non-singleton dimension.
%
%    S = MEDIAN(X,DIM) takes the median along the dimension DIM. 
%
%    S = MEDIAN(X, 0) does the same as above, but returns the total median per object.
%
% Example: s=estruct(-10:10); median(s,0) == 0
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/uminus, estruct/abs, estruct/real, estruct/imag, estruct/uplus

b = unary(a, mfilename, varargin{:});

if nargin > 1 && isequal(varargin{1},0)
  if iscell(b)
    b = cell2mat(b);
  end
end
