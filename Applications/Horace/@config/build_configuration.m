function build_configuration(this, default_config_fun, config_name)
% Construct a configuration structure, save in memory and to file.
%
%   >> build_configuration (config, default_config_fun, config_name)
%
% *** NOTE: This method has to be public because it is used by the constructor of all
%           configuration objects other than the root configuration object. It should
%           not be used in any other context.
%
% Input:
% ------
%   config              Instance of the root configuration class
%   default_config_fun  Function that returns default structure for the configuration class
%   config_name         Name of the configuration class
%
%
% In detail:
% ----------
% Set a configuration in memory from previously saved configuration structure stored
% in a .mat file associated with the particular configuration class name. 
%
% If the file does not exist (e.g. not been created before or has been deleted), or the
% contents of the file are out of date (as determined from the default structure
% constructor), then the default structure is saved to file.
%
% In either case, the current configuration and default configuration are saved in
% memory.

% $Revision: 261 $ ($Date: 2013-09-22 14:56:19 +0100 (Sun, 22 Sep 2013) $)


% Get root configuration name
root_config_name=mfilename('class');

% Forbid this function being called on the root configuration class
if strcmpi(root_config_name,config_name)
    error('Cannot call build_configuration on the root configuration object')
end

% Determine if the class has already been constructed - if so, then nothing to do
if config_store(config_name)
    return
end

% Get default configuration - and make some checks that there are no developer errors
default_config_data=default_config_fun();       % () is required to indicate that this is a function, even though it takes no arguments
if isstruct(default_config_data)
    [valid,mess,default_config_data]=check_fields_valid(default_config_data,root_config_name);  % To check no developer errors
    if ~valid, error('Fields not all valid in default configuration: %s',mess), end
else
    error('Default configuration is not a structure - developer error')
end

% Get stored configuration, if any
file_name = config_file_name (config_name);
[saved_config_data,ok,mess] = load_config (file_name);

% Build configuration from file, if can. Note that load_config will return
% ok==true if the config file does not exist but with config_data_saved empty.
% We check the validity of the fields in the config file because it may have
% been constructed via a route other than the constructor for this configuration.
if ~isempty(saved_config_data)   % configuration data read from file
    % Check fields in saved configuration match those in child default
    if isequal(fieldnames(saved_config_data),fieldnames(default_config_data))
        [valid,mess,saved_config_data_out]=check_fields_valid(saved_config_data,root_config_name);
        if valid
            if ~isequal(saved_config_data_out,saved_config_data)
                warning('Fields are valid but will be updated and saved to match new format')
                [ok,mess]=save_config(file_name,saved_config_data_out);
                if ~ok, error(mess), end
            end
            config_store(config_name,saved_config_data_out,default_config_data,this)
            return
        else
            warning(['Fields not all valid in saved configuration for %s: %s.',...
                '\nIt will be updated with default values.'],config_name,mess)
        end
    else
        warning('CONFIG:build_configuration','Out of date configuration format for %s.\nIt will be updated with default values.',config_name)
    end
elseif ~ok
    warning('CONFIG:build_configuration','%s \n Building configuration for %s with default values',mess,config_name)
end

% Save configuration from defaults.
[ok,mess]=save_config(file_name,default_config_data);
if ~ok, error(mess), end
config_store(config_name,default_config_data,default_config_data,this);
