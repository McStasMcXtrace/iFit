function S=get_all(this,opt)
% Retrieve all initialised configurations
%
%   >> S = get_all (config_obj)             % retrive all current configurations
%   >> S = get_all (config_obj,'defaults')  % retrive all default configurations
%
% Input
% -----
%   config_obj  Configuration object for which all configuration objects
%               that share the same root configuration object will be found, or 
%               the root configuration object itself.
%   
% Keyword:
%   'defaults'  If present, then retrieve the default configuration values for those
%              conmfiguration objects that are initialised.
%               If not present, then retrieve the current values
%
% Output:
% -------
%   S           Structure with fields matching names of all configurations that share the
%               same root configuration object as the input object, and each
%               field contains a structure whose fields are those of the configuration.

% $Revision: 261 $ ($Date: 2013-09-22 14:56:19 +0100 (Sun, 22 Sep 2013) $)

if nargin==1
    S=config_store(false);
elseif nargin==2 && ischar(opt) && ~isempty(opt) && size(opt,1)==1 && strncmpi(opt,'defaults',length(opt))
    S=config_store(true);
else
    error('Check input arguments')
end
