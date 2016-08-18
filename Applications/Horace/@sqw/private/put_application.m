% Write application information data strcuture to file
%
%   >> mess = put_application (application)
%
% Input:
%   fid             File identifier of output file (opened for binary writing)
%   application     Data structure with fields below
%
% Output:
%   mess            Message if there was a problem writing; otherwise mess=''
%
%
% Fields written to file are:
%   application.name        Name of application that wrote the file
%   application.version     Version number of the application
%