% Get indicies that define ranges of contiguous elements from an n-dimensional
% array of bins of elements, where the bins partially or wholly lie
% inside a hypercuboid volume that on the first three axes can be rotated and
% translated w.r.t. to the hypercuboid that is split into bins.
%
% Output will not necessarily be strictly contiguous blocks, as the routine handles the first
% three dimensions separately from the following ones. The blocks are contiguous within
% the first three dimensions, however.
%
%   >> [nstart,nend] = get_nrange (urange,rot,trans,nelmts,,p1,p2,p3,...)
% 
% Input:
%   urange(2,n)     Range to cover: 2 x n array of [urange_lo; urange_hi]
%   rot             Matrix [3x3]     --|  that relate a vector expressed in the
%   trans           Translation [3x1]--|  frame of the bin boundaries to those of urange:
%                                             r'(i) = A(i,j)(r(j) - trans(j))
%                                  (trans is the vector from the origin of the frame
%                                   in which the bins are expressed to that in which
%                                   urange is expressed)
%   nelmts          Array of number of points in n-dimensional array of bins
%                       e.g. 3x5x7 array such that nelmts(i,j,k) gives no. points in (i,j,k)th bin
%   p1(:)           Bin boundaries along first axis
%   p2(:)           Similarly axis 2
%   p3(:)           Similarly axis 3
%    :                      :
%
% Output:
%   nstart          Array of starting values of contiguous block in nelmts(:)
%   nend            Array of finishing values
%
%   If the region defined by urange lies outside the bins, or there are no elements in the range (i.e. the
%   bins that are in the range contain no elements) then nstart and nend returned as empty array [].
%