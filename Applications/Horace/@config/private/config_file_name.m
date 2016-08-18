function file_name = config_file_name (config_name)
% Name of file containing configuration data
%
%   >> file_name = config_file_name (config_name)
%
% Input:
% -----
%   config_name     Name of configuration object
%
% Output:
% -------
%   file_name       Name of file containing the stored value of the
%                   named configuration

% $Revision: 261 $ ($Date: 2013-09-22 14:56:19 +0100 (Sun, 22 Sep 2013) $)

%--> The block to provide compatibility between matlab 2008a and 2007b where
% mfilename behaviour changes
[fd,ff]=fileparts(mfilename('class'));
if isempty(fd) 
    root_config_name = ff;
else
    root_config_name = fd;
end
%<--
fetch_default=false;
config_data = config_store(root_config_name,fetch_default);
file_name=fullfile(config_data.config_folder_path, [config_name,'.mat']);
