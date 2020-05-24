function a = flipud(a)
%  FLIPUD Flip object Signal in up/down direction.
%    FLIPUD(X) returns X with columns preserved and rows flipped
%    in the up/down direction.
%    With 2D data sets, the axis rank 1 is inverted.
%
% Example: s=iData(rand(3)); all(all(flipud(s)==flipud(s{0})))
% Version: $Date$ $Version$ $Author$
% See also iData, iData/fliplr, iData/flipdim

a = unary(a, 'flipdim',1);

