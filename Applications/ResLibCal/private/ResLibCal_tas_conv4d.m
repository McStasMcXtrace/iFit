function signal = ResLibCal_tas_conv4d(data, config, frame)
% model=ResLibCal_tas_conv4d(data or model, tas_config, frame) 4D convolution 
%   function for  neutron Triple-Axis Spectrometers.
%
% This function builds a fit model from:
%   * a dispersion model (iFunc)
%     the dispersion should follow the syntax: dispersion(p, H,K,L,E) where
%     p holds the dispersion model parameters, and H,K,L,E are the 4D coordinates
%     where to evaluate the dispersion (one single mode).
% OR* a data set (iData) which contains some of ResCal parameters, Model, ...
%
%   * a TAS configuration, given as a ResLib/EXP structure, or a configuration
%       file (ResLibCal, ResCal, ResTrax, ResCal5, aFIT/SPEC). If config is left 
%       as empty ('') or not given,  the current ResLibCal configuration will be used.
%
%   * the reference frame to use for the dispersion, as the 'spec' spectrometer 
%     one (x //Q) or the lattice one 'rlu' with a* defined for the lattice orientation.
%
%   IMPORTANT: the H,K,L axes are those defined by the reference frame.
%   When the frame is set as the spectrometer frame 'spec', the 'H' direction is 
%     the Q longitudinal, 'L' is vertical, and 'K' is transverse.
%   When the frame is set as the lattice frame 'rlu', the 'H' direction is along
%     the axis 'a*', the 'K' axis is along b*, and 'H' is along c*.
%   When the frame is set to 'ABC', the 'H' direction is along the user defined 
%     vector 'A' in the scattering plane, the 'K' axis is perpendicular to 'A'
%     in the 'AB' scattering plane, and 'L' is vectical.
%
% The model can then be used for refinement as a usual fit model:
%   fits(model, data_set, parameters, options, constraints)
%
% Example:
% s=sqw_vaks('KTaO3');    % create a 4D S(q,w) perovskite model
% t=ResLibCal(s);        % convolute it with a TAS resolution, and open ResLibCal.
% w=linspace(0.01,20,50); qh=0.6*ones(size(w)); qk=0*qh; ql=qk; 
% signal =iData(t, [], qh,qk,ql,w);
% signal2=iData(s, [], qh,qk,ql,w);
% figure; plot(squeeze([signal signal2*100])); % plot the dispersion and simulated measurement
% % now plot the 4D dispersion with the scan in it
% qh=linspace(0.01,.7,50);qk=qh; ql=qh; w=linspace(0.01,10,51);
% f=iData(s,[],qh,qk,ql,w);
% figure; surf(log(f(:,:,1,:)),'half'); hold on; scatter3(log(signal(:,:,1,:)),'filled');

signal = [];

if nargin < 1
  data = '';
end
if nargin < 2
  config = '';
end
if nargin < 3
  frame = '';
end

if isempty(data)
  error([ mfilename ': no dispersion/data set given for convolution.' ]);
end

% CHECKS and information before convolution ====================================

% open ResLibCal;
ResLibCal('silent');

% configuration is empty, or ResLibcal: will use current config
if ischar(config) && isdir(config), config=''; end

if isempty(config) || (ischar(config) && any(strcmp(lower(config),{'tas','rescal','reslib','reslibcal'})))
  % get the current configuration
  config = ResLibCal('silent','compute');
elseif isstruct(config) && (isfield(config,'EXP') || isfield(config, 'method') || isfield(config, 'ResCal'))
  % load ResLibCal (structure) configuration
  try
    config = ResLibCal('silent',config);
  end
elseif ~isempty(dir(config)) && ~isdir(config)
  % the configuration is a file: we get the configuration and load it into ResLibCal GUI.
  config=ResLibCal('silent',config);
end

% handle array of models/data sets
if numel(data) > 1
  signal = [];
  for index=1:numel(data)
    signal = [ signal ResLibCal_tas_conv4d(data(index), config, frame) ];
  end
  return
end

if isempty(frame), 
  if ndims(data) == 4, frame='rlu'; 
  else frame='spec'; end
end



% assemble the convoluted model/data set
if isa(data,'iFunc')
  % display message about [abc] or [xyz] cloud usage, and axes used
  disp([ mfilename ': Using reference frame: ' frame ])
  if isstruct(config) && isfield(config, 'resolution')
    if ~iscell(config.resolution), resolution={ config.resolution }; 
    else resolution = config.resolution; end
    disp(resolution{1}.(frame).README);
  end
  signal = ResLibCal_tas_conv4d_model(data, config, frame);
elseif isa(data,'iData')
  signal = ResLibCal_tas_conv4d_data(data, config, frame);
end

% ==============================================================================

function signal = ResLibCal_tas_conv4d_model(dispersion, config, frame)
% build the conv(model) using ResLibCal resolution 'cloud'
%
% the model will be used as:
%    model(p, H,K,L,W)
% or
%    model(p, scan_axis)
% or
%    model(p, iData_with_scan_axis)

% select computation mode: normal or fast
% the 'hkle' mode uses the ResLibCal cloud as is, i.e. HKLE with NMC points (e.g. 200)
% the 'hkl'  mode only uses the HKL cloud projection, and builds a E regular axis
computation_mode = 'hkl';

if isempty(dispersion), signal = []; return; end

% assemble the convoluted model ================================================

if ndims(dispersion) ~= 2 && ndims(dispersion) ~= 4
  disp([ mfilename ': the Model should be 2D or 4D, not ' num2str(ndims(dispersion)) '. Skipping.' ])
  signal = dispersion;
  return
end

% check if the model is not already convoluted with TAS
if any(~cellfun(@isempty,strfind(dispersion.Expression,'ResLibCal')))
  signal = dispersion;
  return
end

% the built model parameters will be that of the model
signal.Parameters = dispersion.Parameters;
signal.ParameterValues = dispersion.ParameterValues;

if ~isempty(config)
  if ~ischar(config)
    configs = class2str(config,'eval');
    if length(configs) > 20, configs=[ configs(1:20) '...' ]; end
  else configs = config; end
  signal.Name       = [ 'conv(' dispersion.Name ', ResLibCal(''' configs ''')) [' mfilename ']' ];
  signal.Description= [ '(' dispersion.Description ') convoluted by (TAS 4D resolution function)' ];
else
  signal.Name       = [ 'conv(' dispersion.Name ', ResLibCal)' ];
  signal.Description= [ '(' dispersion.Description ') convoluted by (TAS 4D resolution function with configuration ' config ')' ];
end
signal.Dimension = dispersion.Dimension;
signal.Guess     = dispersion.Guess;
signal.Constraint= dispersion.Constraint;

% we store the dispersion into UserData so that we can evaluate it at feval
signal.UserData.dispersion = dispersion;
signal.UserData.config     = config;

% create the Expression...

% TODO: integration can be achieved using a Gauss-Hermite estimate, which is much
% faster, as in B. Hennion AFIT.

% if dispersion.Dimension == 2 (liquid,powder,gas,glass,polymer...), then use |q|,w as axes. Should use spec frame
if dispersion.Dimension == 2
  % TODO: must use reciprocal space frame, as in sqw_powder
  plugin = 'cloud{1} = sqrt(cloud{1}.^2+cloud{2}.^2+cloud{3}.^2); cloud(2:3)=[]; ';
else plugin='';
end

% handle normal and fast computation mode.
if ~strfind(computation_mode, 'hkle')
  % transform energy random sampling into a regular one
  plugin = [ plugin 'cloud = { cloud{1:3} linspace(min(cloud{4}),max(cloud{4}),ceil(numel(cloud{1})^.33)) };' ];
end

% we use the Sqw evaluation as:
% kpath 2D: [ x y z ]     as vectors, same length, same orientation, [ t ] not same length or orientation
signal.Expression = { ...
  '% check if config is given as additional argument during evaluation of model', ...
  'if numel(varargin) >= 1 && ~isempty(varargin{1}) && (isstruct(varargin{1}) || ~isempty(dir(varargin{1})))', ...
  '  config = varargin{1};', ...
  'else config=[]; end', ...
  'if isempty(config) && ~isempty(this.UserData.config), config = this.UserData.config; end', ...
  '% save current HKLE location', ...
  'hkle = ResLibCal(''silent'',''hkle'');', ...
  '% put it in a try catch so that we always restore the HKLE location', ...
  'try', ...
  '% compute the resolution, using either given config, or default/GUI', ...
  'if isempty(config), out=ResLibCal(''silent'', x,y,z,t);', ...
  'else out=ResLibCal(''silent'', config, x,y,z,t); end', ...
  '% get the MC cloud', ...
  'if ~iscell(out.resolution), resolution={ out.resolution };', ...
  'else resolution = out.resolution; end', ...
  'signal=zeros(size(resolution));', ...
  'dispersion=this.UserData.dispersion;', ...
  'for index=1:numel(resolution)', ...
  '  if ~resolution{index}.R0, continue; end', ...
[ '  cloud=resolution{index}.' frame '.cloud;' ], ...
     plugin, ...
  '  [this_signal,dispersion]=feval(dispersion, p, cloud{:});', ...
  '  dispersion.ParameterValues = p;', ...
  '  signal(index) = sum(this_signal(:))*resolution{index}.R0/numel(cloud{1});', ...
  'end % for', ...
  'this.UserData.dispersion = dispersion;', ...
  'end % try', ...
  '% restore previous HKLE location', ...
  'if ~isempty(hkle) ResLibCal(''silent'',hkle{:}); end'};
  
signal = iFunc(signal);
disp(signal.Name);

% ==============================================================================

function c = ResLibCal_tas_conv4d_data(a, config, frame)

% convolute the iData.Model with ResLibCal, overlay parameters from the Data set
% and provide missing axes from the data set into the Model

% search for missing axes (e.g. 1D -> 4D)
% we must create a new 4D object'a' which has proper axes
axes_symbols = {'QH','QK','QL','EN'};
sz = size(a.Signal);
for index = 1:numel(axes_symbols) % also searches for 'lower' names
  % the axis must be either scalar, vector with size(Signal, rank), or same size as the Signal
  options = 'exact numeric';
  if index > 1, options = [ options ' cache' ]; end
  [match, types, nelements]=findfield(a, axes_symbols{index}, options);
  % must get the longest field (prefer column from data file rather than simple
  % static scalar)
  if ischar(match) match= { match }; end
  if ~isempty(match)
    % first sort matches with the number of elements
    [nelements, index_sort] = sort(nelements);
    index_sort = index_sort(end:-1:1);  % sort by decreasing size
    match = match(index_sort);
    % look for axes symbols by 
    for index_match = 1:numel(match)
      this_axis = get(a, match{index_match});
      if isscalar(this_axis) ...
        || (numel(sz) >= index && numel(this_axis) == sz(index)) ...
        || (numel(this_axis) == prod(sz))
        a=setaxis(a, index, match{index_match});
        a=label(a, index, axes_symbols{index});
        disp([ mfilename ': setting axis ' num2str(index) ' ' axes_symbols{index} ' as ' match{index_match} ]);
        break
      end
    end
  end
end

% *** RESCAL parameters ***
if isfield(config,'ResCal')
  rescal = config.ResCal; % from the interface
else
  rescal = [];
end

% search for ResCal parameters in the data set
if ~isempty(rescal) && isstruct(rescal)
  % make sure the HKLE coordinates are single values
  for index = 1:numel(axes_symbols)
    if isfield(rescal, axes_symbols{index})
      rescal.(axes_symbols{index}) = mean(rescal.(axes_symbols{index}));
    end
  end
  disp([ mfilename ': using current parameters from the main ResLibCal window.' ])
  disp('    To upload parameters from the data file:')
  disp([ '      ' char(a) ])
  disp('    open it with File>Open menu item ResLibCal');
end

ResLibCal('silent','open',rescal);


% get the model from the data set
model = [];
if isfield(a, 'Model')
  model = get(a, 'Model');
elseif ~isempty(findfield(a, 'Model'))
  model = get(a, findfield(a, 'Model', 'cache first'));
end

% update the embedded model with 4D convolution
if isa(model, 'iFunc') && ~isempty(model) && any(ndims(model) == [2 4])

  % 4D convolution/ResLibCal
  model = ResLibCal(model);  % calls iFunc.conv == ResLibCal(model)
  % the ResLibCal config is stored in model.UserData.config
  
  % store the ResCal parameters so that they can be used in subsequent convolutions
  % the parameters are automatically used in the model Expression 
  % which calls "ResLibCal('silent', config, x,y,z,t)" in model.Expression
  % where config = this.UserData.config which contains ResCal
  for f=fieldnames(rescal)'
    model.UserData.config.ResCal.(f{1}) = rescal.(f{1});
  end
  
  % to force the use of Rescal parameters, and not ResLib ones, we remove the 
  % 'EXP' member from ResLibCal config.
  % update the Model.UserData stuff with this so that ResLibCal can use it later
  if isfield(model.UserData.config, 'EXP')
    model.UserData.config = rmfield(model.UserData.config, 'EXP');
  end
  
  % update new model
  a=setalias(a, 'Model',model);
end

c = copyobj(a); % finally make a new object
