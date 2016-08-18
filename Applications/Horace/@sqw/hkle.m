% [h,k,l,e] for points in the coordinates of the display axes for an sqw object from a single spe file
%
%   >> [qe1,qe2] = hkle(w,x)
%
% Input:
% ------
%   w       sqw object
%   x       Vector of coordinates in the display axes of an sqw object
%           The number of coordinates must match the dimensionality of the object.
%          e.g. for a 2D sqw object, x=[x1,x2], where x1, x2 are column vectors.
%           More than one point can be provided by giving more rows
%          e.g.  [1.2,4.3; 1.1,5.4; 1.32, 6.7] for 3 points from a 2D object.
%           Generally, an (n x nd) array, where n is the number of points, nd
%          the dimensionality of the object.
%
% Output:
% -------
%   qe1     Components of momentum (in rlu) and energy for each bin in the dataset
%           Generally, will be (n x 4) array, where n is the number of points
%
%   qe2     For the second root
%%   Overloaded methods:
%      sqw/hkle
%      sqw/hkle
%