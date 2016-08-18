% Read a collection of sqw files and sort the pixels from those files onto a common grid.
% Write the results to disk, one file per input sqw file.
%
%   >> [grid_size, urange] = write_nsqw_to_nsqw (dummy, infiles, outfiles)
%   >> [grid_size, urange] = write_nsqw_to_nsqw (dummy, infiles, outfiles, grid_size_in)
%
% Input:
%   dummy           Dummy sqw object  - used only to ensure that this service routine was called
%   infiles         Cell array or character array of file name(s) of input file(s)
%   outfiles        Cell array or character array of full name(s) of output file(s)
%   grid_size_in    [Optional] Scalar or row vector of grid dimensions.
%                  Default is [10,10,10,10]
%   urange_in       [Optional] Range of data grid for output. If not given, then uses smallest hypercuboid
%                  that encloses the whole data range.
%
% Ouput:
%   grid_size       Actual grid size used (size is unity along dimensions
%                  where there is zero range of the data points)
%   urange          Actual range of grid
%
%%   Overloaded methods:
%      sqw/write_nsqw_to_nsqw
%      sqw/write_nsqw_to_nsqw
%