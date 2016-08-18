% Get indicies that define ranges of contiguous bins that partially or wholly lie
% inside a cuboid volume that is rotated and translated w.r.t. to the cuboid that is split
% into bins.
% Works for non-uniformly spaced bin boundaries.
%
%   >> [istart,iend] = get_irange_rot(urange,p1,p2,p3)
%
% Input:
%   urange(2,3)     Range to cover: 2 x 3 array of [urange_lo; urange_hi]
%                  in units of a rotated and shifted coordinate frame
%   rot             Matrix     --|  that relate a vector expressed in the
%   trans           Translation--|  frame of the bin boundaries to those of urange:
%                                      r'(i) = A(i,j)(r(j) - trans(j))
%                                  (trans is the vector from the origin of the frame
%                                   in which the bins are expressed to that in which
%                                   urange is expressed)
%   p1(:)           Bin boundaries along first axis
%   p2(:)           Similarly axis 2
%   p3(:)           Similarly axis 3
% 
% Ouptut:
%   istart          Bin index values. If range is outside the bins then returned as empty
%   iend            Upper index values. If range is outside the bins then returned as empty
%
% Note:
%   The algorithm is an n^3 algorithm - good for small grids, but could be improved for
%   for large grids.
%