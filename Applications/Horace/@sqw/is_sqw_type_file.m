% Determine if file contains data of an sqw-type object or dnd-type sqw object
% 
%   >> [sqw_type, ndims, data_source, mess] = is_sqw_type_file(sqw, infile)
%
% Input:
% --------
%   w           Dummy sqw object, whose sole purpose is to get to this function
%   infile      File name, character array of file names, or cellstr of file names
%
% Output:
% --------
%   sqw_type    =true  if sqw-type contents; =false if dnd-type contents (array)
%   ndims       Number of dimensions (array if more than one file)
%   filename    Cell array of file names (even if only one file, this is still a cell array)
%   mess        Error message; blank if no errors, non-blank otherwise
%%   Overloaded methods:
%      sqw/is_sqw_type_file
%      sqw/is_sqw_type_file
%