function a = flipud(a, varargin)
%  FLIPDIM Flip matrix along specified dimension.
%    FLIPDIM(X,DIM) returns X with dimension DIM flipped.
%    FLIPDIM(X,1) is equivalent to FLIPUD
%    FLIPDIM(X,2) is equivalent to FLIPLR
%
% Example: s=iData(rand(3,4,5)); z=all(flipdim(s,3)==flipdim(s{0},3)); unique(z)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/fliplr, fliplr, iData/flipud, flipud

a = unary(a, 'flipdim', varargin{:});

