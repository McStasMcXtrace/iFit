% the function verifies the grid sizes, sets defaults if grid size is wrong
% and fills in the bins boundaries for grid axis
%
% inputs:
% grid_size_in  --  initial grdis size guess
% urange        --  range of the input data
% nd            --  number of data dimentions; should be equal to
%                   size(grid_size,2) if grid_size is defined. defines the
%                   grid_size otherwise;
% outputs:
% grid_size     -- verified or constructed grid size
% p             -- bin coordinates;
%
% Construct grid_size array if necessary
%
% $Revision: 601 $ ($Date: 2012-02-08 14:46:10 +0000 (Wed, 08 Feb 2012) $)
%
%