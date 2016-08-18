% Write the main header block for the results of performing calculate projections on spe file(s).
%
%   >> mess = put_sqw_header (fid, data)
%
% Input:
%   fid             File identifier of output file (opened for binary writing)
%   data            Data structure which must contain (at least) the fields listed below
%
% Output:
%   mess            Message if there was a problem writing; otherwise mess=''
%
%
% Fields written to file are:
%   data.filename   Name of file that was the source of sqw data structure, excluding path
%   data.filepath   Path to file including terminating file separator
%   data.title      Title of sqw data structure
%   data.nfiles     Number of spe files that contribute
%
%
% Notes:
% ------
%   There are some other items written to the file to help when reading the file using get_sqw_data. 
% These are indicated by comments in the code.
%