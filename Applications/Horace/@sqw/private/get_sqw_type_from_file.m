% Get sqw_type and dimensionality of an sqw file on disk
%
% Syntax:
%   >> [sqw_type, ndims, mess] = get_sqw_type_from_file(infile)
%
% Input:
% --------
%   infile      File name
%
% Output:
% --------
%   sqw_type    =true  if sqw-type contents; =false if dnd-type contents
%   ndims       Number of dimensions
%   mess        Error message; blank if no errors, non-blank otherwise
%