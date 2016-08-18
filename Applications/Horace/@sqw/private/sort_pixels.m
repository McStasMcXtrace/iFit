% Reorder the pixels according to increasing bin index in a Cartesian grid.
%
%   >> [ix,npix,p,grid_size]=sort_pixels(u,urange,grid_size_in)
%
% (In the following, nd=no. dimensions, npix_in=no. pixels on input, npix=no. pixels in urange)
% Input:
%   u               [nd x npix_in] array of the coordinates in Cartesian grid
%   urange          Range for the grid (2 x nd)
%   grid_size_in    Scalar or row vector (1 x nd) of number of bins along each axis
%
% Output: 
%   ix(npix,1)      Index array by which the pixels have been reordered
%                  i.e. u1(ix) gives the reordered values of u1
%   npix(nbin,1)    Number of contributing pixels to the bins in the Cartesian grid
%                  as a column vector. Bin indicies of reshaped Cartesian grid
%                  are numbered in the same sequence as would be returned by the
%                  matlab instrinsic sub2ind)
%   p               Cell array [1 x nd] of column vectors of bin boundaries
%                  in Cartesian grid
%   grid_size       Scalar or row vector (1xnd) of number of actual bins along each axis
%                  This may differ from the input grid_size if the range along any of the
%                  axes is zero: in this case the size of the grid along those axes = 1
%   ibin(npix,1)    Column vector with list of bins to which the sorted pixels contribute
%                  Available for convenience as it can be constructed from npix:
%                       ibin(1:nbin(1))=1
%                       ibin(nbin(1):nbin(2))=2
%                               :
%                       ibin(nbin(end-1):nbin(end))=length(nbin)
%
%