function y = tas_conv4d(dispersion, config, frame) 
% model=tas_conv4d(dispersion_model, tas_config) 4D convolution function for 
%   neutron Triple-Axis Spectrometers.
%
% This function buildsa fit model from:
%   * a dispersion model (iFunc or expression)
%     the dispersion should follow the syntax: dispersion(p, H,K,L,E) where
%     p holds the dispersion model parameters, and H,K,L,E are the 4D coordinates
%     where to evaluate the dispersion (one single mode).
%
%   * a TAS configuration, given as a ResLib/EXP structure, or a configuration
%       file (ResLibCal, ResCal, ResTrax, ResCal5, aFIT/SPEC). If config is left 
%       as empty ('') or not given,  the current ResLibCal configuration will be used.
%
%   * the reference frame to use for the dispersion, as the 'xyz' spectrometer 
%     one (x //Q) or the lattice one 'abc' with A defined for the lattice orientation.
%
%   IMPORTANT: the H,K,L axes are those defined by the reference frame.
%   When the frame is set as the spectrometer frame 'xyz', the 'H' direction is 
%     the Q longitudinal, 'L' is vertical, and 'K' closes the direct frame.
%   When the frame is set as the lattice frame 'abc', the 'H' direction is along
%     the axis 'A' defined for the lattice orientation, the 'K' axis is orthogonal
%     in plane, and 'L' is vertical.
%
%
% The model can then be used for refinement as a usual fit model:
%   fits(model, data_set, parameters, options, constraints)

if nargin < 1
  dispersion = '';
end
if nargin < 2
  config = '';
end
if nargin < 3
  frame = '';
end
if isempty(frame), frame='abc'; end
frame = lower(frame);

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

% check if config is given as additional argument during evaluation of model'
'if numel(varargin) >= 1 && isstruct(varargin{1})'
'  config = varargin{1};'
'else'
if ~isempty(config) && isstruct(config)
if isfield(config,'EXP'), EXP=config.EXP; else EXP=config; end
disp('Using given static configuration as default.')
[ class2str('config', EXP) ]; % store the ResLib config
else
disp('Using current ResLibCal configuration as default.')
'config='''';';
end
'end'

% display message about [abc] or [xyz] cloud usage, and axes used ?
disp([ 'Using reference frame: ' frame ])
if isstruct(config) && isfield(config, 'resolution')
  if ~iscell(config.resolution), resolution={ config.resolution }; 
  else resolution = config.resolution; end
  end
  disp('Frame axes:')
  resolution{1}.(frame).FrameStr
end

% when starting, auto-update is set to off to avoid long computation+display
ResLibCal('autoupdate','off');

% store the dispersion
% it may have a zone-center HKL location
% the built model parameters will be that of the model

% compute the resolution, using either given config, or default/GUI
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
'signal=zeros(size(resolution));'
'for index=1:numel(resolution)'
'cloud=resolution{index}.' frame '.cloud;'
% if dispersion.Dimension == 2 (liquid,powder,gas,glass,polymer...), then use |q|,w as axes
if dispersion.Dimension == 2
'cloud{1} = sqrt(cloud{1}.^2+cloud{2}.^2+cloud{3}.^2); cloud(2:3)=[];'
end
'signal(end+1)=sum(feval(dispersion, p, cloud{:}))*resolution{index}.R0/numel(cloud{1})'
'end'


