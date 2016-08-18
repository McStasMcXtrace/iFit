% Load header information of VMS format ASCII .spe file
%   >> data = get_spe(filename)
%
% data has following fields:
%   data.filename   Name of file excluding path
%   data.filepath   Path to file including terminating file separator
%   data.ndet       Number of detector groups
%   data.en         Column vector of energy bin boundaries
%