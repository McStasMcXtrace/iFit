function b = mean(a, varargin)
%  MEAN   Average or mean value of Signal.
%    S = MEAN(X) is the mean of the Signal of the object X. If
%    X is a N-D object, mean(X) operates along the first
%    non-singleton dimension.
%
%    S = MEAN(X,DIM) takes the mean along the dimension DIM. 
%
%    S = MEAN(X, 0) does the same as above, but returns the total mean per object.
%
% Example: s=iData(-10:10); mean(s,0) == 0
% Version: $Date$ $Version$ $Author$
% See also iData, iData/uminus, iData/abs, iData/real, iData/imag, iData/uplus

b = unary(a, mfilename, varargin{:});

if nargin > 1 && isequal(varargin{1},0)
  if iscell(b)
    b = cell2mat(b);
  end
end
