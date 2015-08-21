function signal = tas_conv4d(dispersion, config, frame)
% model=tas_conv4d(dispersion_model, tas_config) 4D convolution function for 
%   neutron Triple-Axis Spectrometers.
%
% This function builds a fit model from:
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
% The model can then be used for refinement as a usual fit model:
%   fits(model, data_set, parameters, options, constraints)
%
% Example:
%   s=sqw_vaks('KTaO3');  % create a 4D S(q,w) perovskite model
%   t=tas_conv(s);        % convolute it with a TAS resolution, and open ResLibCal.

if nargin < 1
  dispersion = '';
end
if nargin < 2
  config = '';
end
if nargin < 3
  frame = '';
end

frame = lower(frame);
if isempty(dispersion)
  error([ mfilename ': no dispersion given for convolution' ]);
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

% CHECKS and information before building the model =============================

% open ResLibCal;
ResLibCal;

% configuration is empty, or ResLibcal: will use current config
if isdir(config), config=''; end

if ischar(config) && any(strcmp(lower(config),{'tas','rescal','reslib','reslibcal'}))
  % get the current configuration
  config = ResLibCal;
elseif isstruct(config) && (isfield(config,'EXP') || isfield(config, 'method'))
  % load ResLibCal (structure) configuration
  try
    config = ResLibCal(config);
  end
elseif ~isempty(dir(config)) && ~isdir(config)
  % the configuration is a file: we get the configuration and load it into ResLibCal GUI.
  config=ResLibCal(config);
end

% when starting, auto-update is set to off to avoid long computation+display
disp([ mfilename ': setting ResLibCal:autoupdate to OFF for efficiency in computations.' ]);
ResLibCal('autoupdate','off');

if isempty(frame), 
  if dispersion.Dimension == 4, frame='abc'; 
  else frame='xyz'; end
end

% display message about [abc] or [xyz] cloud usage, and axes used
disp([ mfilename ': Using reference frame: ' frame ])
if isstruct(config) && isfield(config, 'resolution')
  if ~iscell(config.resolution), resolution={ config.resolution }; 
  else resolution = config.resolution; end
  disp('Frame axes:')
  resolution{1}.(frame).FrameStr
end

% assemble the convoluted model ================================================

% the built model parameters will be that of the model
signal.Parameters = dispersion.Parameters;
if ~isempty(config)
  signal.Name       = [ 'conv(' dispersion.Name ', ResLibCal(''' config ''')) [' mfilename ']' ];
  signal.Description= [ '(' dispersion.Description ') convoluted by (TAS 4D resolution function)' ];
else
  signal.Name       = [ 'conv(' dispersion.Name ', ResLibCal)' ];
    signal.Description= [ '(' dispersion.Description ') convoluted by (TAS 4D resolution function with configuration ' config ')' ];
end
signal.Dimension = dispersion.Dimension;
signal.Guess     = dispersion.Guess;

% we store the dispersion into UserData so that we can evaluate it at feval
signal.UserData.dispersion = dispersion;
signal.UserData.config     = config;

% create the Expression...
% TODO: in a fit procedure, as the coordinates xyzt (HKLE) do not change, the clouds 
% may only be computed once. If cached, the subsequent fit steps will be faster.
% then the function should be able to send its clouds out side, e.g. p='cloud'.
% Then the fit function can send this cache as additional argument during the
% fit.
% TODO: integration can be achieved using a Gauss-Hermite estimate, which is much
% faster

% if dispersion.Dimension == 2 (liquid,powder,gas,glass,polymer...), then use |q|,w as axes. Should use xyz frame
if dispersion.Dimension == 2
  liq = 'cloud{1} = sqrt(cloud{1}.^2+cloud{2}.^2+cloud{3}.^2); cloud(2:3)=[];';
else liq='';
end
signal.Expression = { ...
  '% check if config is given as additional argument during evaluation of model', ...
  'if numel(varargin) >= 1 && (isstruct(varargin{1}) || ~isempty(dir(varargin{1})))', ...
  '  config = varargin{1};', ...
  'else config=[]; end', ...
  '% compute the resolution, using either given config, or default/GUI', ...
  'if isempty(config), out=ResLibCal(x,y,z,t);', ...
  'else out=ResLibCal(config, x,y,z,t); end', ...
  '% get the [abc] reference cloud', ...
  'if ~iscell(out.resolution), resolution={ out.resolution };', ...
  'else resolution = out.resolution; end', ...
  'signal=zeros(size(resolution));', ...
  'for index=1:numel(resolution)', ...
  '  if ~resolution{index}.R0, continue; end', ...
[ '  cloud=resolution{index}.' frame '.cloud;' ], ...
liq, ...
  '  dispersion=this.UserData.dispersion;', ...
  '  this_signal=feval(dispersion, p, cloud{:});', ...
  '  signal(index) = sum(this_signal(:))*resolution{index}.R0/numel(cloud{1})', ...
  'end' };

signal = iFunc(signal);
