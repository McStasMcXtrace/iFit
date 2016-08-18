% Given an array containing number of points in bins, and a section of
% that array, return column vectors of start and end indicies of contiguous points
% in the column representation of the points.
% Works for any dimensionality 1,2,...
%
%   >> [nstart,nend] = get_nrange(nelmts,irange)
%
% Input:
%   nelmts      Array of number of points in n-dimensional array of bins
%              e.g. 3x5x7 array such that nelmts(i,j,k) gives no. points in (i,j,k)th bin
%   irange      Ranges of section [irange_lo;irange_hi] (assumes irange_lo<=irange_hi)
%              e.g. [1,2,6;3,4,7] means bins 1:3, 2:4, 6:7 along the three axes
% Output:
%   nstart      Array of starting values of contiguous block in nelmts(:)
%   nend        Array of finishing values
%               nstart=nend=[] if there are no elements
%