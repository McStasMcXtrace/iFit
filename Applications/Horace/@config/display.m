function display (this)
% Display a configuration to the screen
%
%   >> display (config_obj)  % config_obj is an instance of a configuration class


% $Revision: 221 $ ($Date: 2012-09-08 09:54:13 +0100 (Sat, 08 Sep 2012) $)

config_data=get(this);
config_name=class(this);
disp(' ')
disp([config_name,' ='])
disp(' ')
disp(config_data)
