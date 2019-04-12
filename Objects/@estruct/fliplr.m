function a = fliplr(a)
%  FLIPLR Flip object Signal in left/right direction.
%    FLIPLR(X) returns X with row preserved and columns flipped
%    in the left/right direction.
%    With 2D data sets, the axis rank 2 is inverted.
%
% Example: s=estruct(rand(3)); all(all(fliplr(s)==fliplr(s{0})))
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/flipud, estruct/flipdim

a = unary(a, 'flipdim',2);

