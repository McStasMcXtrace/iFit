% Determine sqw type from the data component of an sqw data structure
%
%   >> type = sqw_type(data)
%
%   data        data componenet of sqw structure
%   type        ='b','b+' 'a' (valid sqw type) or 'a-' (sqw without pix)
%               ='h' if header part of data structure only
%
%   Simple routine - the assumption is that the data corresponds to a valid type
%   Needs full data structure, not just the header fields
%