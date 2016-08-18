% Calculate qh,qk,ql,en for set of points in an n-dimensional sqw dataset
%
%   >> qw=calculate_qw_points(win)          % (Q,w) when energy is the implicit variable for a direct geometry cut
%   >> [qw1,qw2]=calculate_qw_points(win)   % in general two roots for other cases
%
% Input:
% ------
%   win     Input sqw object created from a single spe file
%
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
%   qw1     Components of momentum (in rlu) and energy for each bin in the dataset
%           Generally, will be (n x 4) array, where n is the number of points
%
%   qw2     For the second root
%
%   If direct geometry, and
%    - energy transfer is the implicit variable to be determined:
%           there is only one root and qw2=[];
%    - if a component of Q is implicit variable to be determined
%           in general there are either zero or two roots
%             - if two roots
%                 qw1 corresponds to the root with more negative component along the infinite integration axis
%                 qw2 corresponds to the more positive component
%             - if no roots
%                 all elements of the corresponding row in qw1, qw2 set to NaN
%
%   If indirect geometry, and
%    - energy transfer is the implicit variable to be determined:
%           can be zero, one or two roots;
%             - if no roots, corresponding row in qw1 and qw2 set to NaN
%             - if there is a root, it will be in qw1, and the row in qw2 set to naN
%             - if two roots, the larger energy trasnfer root is always in qw1
%           
%    - if a component of Q is implicit variable to be determined
%           in general there are either zero or two roots, just as for direct geometry.
%
%
%