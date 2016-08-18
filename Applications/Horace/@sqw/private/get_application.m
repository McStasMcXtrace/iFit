% Read the application block that gives information about the applciation that wrote
% the file
%
% Syntax:
%   >> [application, mess] = get_application (fid, application_in)
%
% Input:
% ------
%   fid             File pointer to (already open) binary file
%   application_in  [optional] Data structure to which the data
%                  fields below will be added or overwrite.
%
% Output:
% -------
%   application     Structure containing fields read from file (details below)
%   mess            Error message; blank if no errors, non-blank otherwise
%
% Fields read from file are:
%   application.name        Name of application that wrote the file
%   application.version     Version number of the application
%