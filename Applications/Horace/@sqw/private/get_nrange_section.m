% Get contiguous ranges of elements from a subsection of an n-dimensional array
%
%   >> [nstart,nend] = get_nrange (urange,nelmts,,p1,p2,p3...)
% 
% Input:
%   urange(2,nd)    Range to cover: 2 x nd array of [urange_lo; urange_hi]
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