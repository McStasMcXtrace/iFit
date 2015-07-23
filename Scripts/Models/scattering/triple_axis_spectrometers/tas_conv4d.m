function y = tas_conv4d(dispersion, config) 
% model=tas_conv4d(dispersion_model, tas_config) 4D convolution function for 
%   neutron Triple-Axis Spectrometers.
%
% This function build a fit model from:
%   * a dispersion model (iFunc or expression)
%     the dispersion should follow the syntax: dispersion(p, H,K,L,E) where
%     p holds the dispersion model parameters, and H,K,L,E are the 4D coordinates
%     where to evaluate the dispersion (one single mode).
%
%   * a TAS configuration, given as a ResLib/EXP structure, or a configuration
%       file (ResLibCal, ResCal, ResTrax, ResCal5, aFIT/SPEC). If config is left 
%       as empty ('') or not given,  the current ResLibCal configuration will be used.
%
% The model can then be used for refinement as a usual fit model:
%   fits(model, data_set, parameters, options, constraints)

if nargin < 1
  dispersion = '';
end
if nargin < 2
  config = '';
end

% config usage:
% if ResLibCal is opened, we send the 'config' to it and prepare for further evaluations.

% if ResLibCal is not opened, we store the configuration in the Expression, and will then
% make calls to out=ResLibCal(config, 'compute');

% ResLibCal must be available (shipped with iFit)

% the model will be used as:
%    model(p, H,K,L,W)
% or
%    model(p, scan_axis)
% or
%    model(p, iData_with_scan_axis)

% the configuration is a file: we get the configuration and load it into ResLibCal GUI.
if ~isempty(dir(config))
  config=ResLibCal(config);
end

% create the Expression...

% store the config, hard coded.
% TODO: the configuration could also be stored as a GLOBAL variable, together with
% the cloud as cache.

% when starting, auto-update is set to off to avoid long computation+display
ResLibCal('autoupdate','off');
if isstruct(config)
[ class2str('config',config) ]
else
'config='''';'
end

% store the dispersion

% compute the resolution
'out=ResLibCal(config, x,y,z,t);'
% TODO: in a fit procedure, as the coordinates xyzt (HKLE) do not change, the clouds 
% may only be computed once. If cached, the subsequent fit steps will be faster.
% the the function should be able to send its clouds out side, e.g. p='cloud'.
% Then the fit function can send this cache as additional argument during the
% fit.
% TODO: integration can be achieved using a Gauss-Hermite estimate, which is much
% faster

% get the [abc] reference cloud
'if ~iscell(out.resolution), resolution={ out.resolution }; else resolution = out.resolution; end'
'signal=[];'
'for index=1:numel(resolution)'
% check if [abc].cloud exists, else use [xyz].cloud
% display message about [abc] or [xyz] cloud usage, and axes used
'cloud=resolution{index}.abc.cloud'
'signal(end+1)=sum(feval(dispersion, p, cloud{:}))*resolution{index}.R0/numel(cloud{1})'
'end'


