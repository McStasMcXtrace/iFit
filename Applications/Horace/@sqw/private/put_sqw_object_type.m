% Write application information data strcuture to file
%
%   >> mess = put_sqw_object_type (sqw_type)
%
% Input:
%   fid             File identifier of output file (opened for binary writing)
%   sqw_type        Type of sqw object: =1 if sqw type; =0 if dnd type
%   ndims           Number of dimensions of sqw object
%
% Output:
%   mess            Message if there was a problem writing; otherwise mess=''
%