% Read the type of sqw object written to file
%
% Syntax:
%   >> [sqw_type, mess] = get_sqw_object_type (fid)
%
% Input:
% ------
%   fid         File pointer to (already open) binary file
%
% Output:
% -------
%   sqw_type    Type of sqw object written to file: =1 if sqw type; =0 if dnd type
%   mess        Error message; blank if no errors, non-blank otherwise
%
%