function a = del2(a, varargin)
%  DEL2 Discrete Laplacian.
%   L = DEL2(U), when U is a matrix, is a discrete approximation of
%   del^2 u = (d^2u/dx^2 + d^2u/dy^2).  The matrix L is the same
%   size as U, with each element equal to the difference between an 
%   element of U and the average of its four neighbors.
%
%   L = DEL2(U), when U is an N-D array, returns an approximation of
%   (del^2 u).
%
%   L = DEL2(U,H), where H is a scalar, uses H as the spacing between
%   points in each direction (H=1 by default).
%
%   L = DEL2(U,HX,HY), when U is 2-D, uses the spacing specified by HX
%   and HY. If HX is a scalar, it gives the spacing between points in
%   the x-direction. If HX is a vector, it must be of length SIZE(U,2)
%   and specifies the x-coordinates of the points.  Similarly, if HY
%   is a scalar, it gives the spacing between points in the
%   y-direction. If HY is a vector, it must be of length SIZE(U,1) and
%   specifies the y-coordinates of the points.
%
%   L = DEL2(U,HX,HY,HZ,...), when U is N-D, uses the spacing given by
%   HX, HY, HZ, etc.
%
% Example: [x,y]=meshgrid(-4:4,-3:3); s=iData(x.*x+y.*y); all(del2(s) == 4)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/gradient, del2, gradient, iData/jacobian

% make sure axes are regularly binned
%a = interp(a);

a = unary(a, 'del2', varargin{:});

