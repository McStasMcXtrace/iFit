function [pars,fval,exitflag,output] = mcstas(instrument, parameters, options, varargin)
% [OPTIMUM,MONITORS,EXITFLAG,OUTPUT] = mcstas(INSTRUMENT, PARAMETERS, OPTIONS, CONSTRAINTS) : run and optimize a McStas simulation
%
% A wrapper to the McStas package to either execute a simulation, or optimize a
%   set of parameters.
%
% input:  INSTRUMENT: name of the instrument description to run (string)
%         PARAMETERS: a structure that gives instrument parameter names and values (structure)
%           values should be given as numerical values, e.g. 2.36
%             use 'constraints' to define min/max and fix some parameters during optimization.
%           values given as strings are used as-is for e.g. string instrument parameter. These can not be optimized.
%         OPTIONS: a structure or string that indicates what to do (structure)
%           options.dir:    directory where to store results (string)
%           options.ncount: number of neutron events per iteration, e.g. 1e5 (double)
%           options.mpi:    number of processors to use with MPI on localhost (integer) 
%           options.seed:   random number seed to use for each iteration (double)
%           options.gravitation: 0 or 1 to set gravitation handling in neutron propagation (boolean)
%           options.compile: 0 or 1 to force re-compilation of the instrument (boolean)
%           options.mode:   'simulate' or 'optimize' (string)
%           options.type:   'minimize' or 'maximize', which is the default (string)
%           options.monitors:  cell string of monitor names, or empty for all (cellstr)
%           options.optimizer: function name of the optimizer to use (string or function handle)
%         as well as other optimizer options such as
%           options.TolFun ='0.1%'   stop when criteria changes are smaller that 0.1%
%           options.Display='final'
%         CONSTRAINTS may be specified as a structure
%           constraints.min=   vector of minimal values for parameters
%           constraints.max=   vector of maximal values for parameters
%           constraints.fixed= vector having 0 where parameters are free, 1 otherwise
%           constraints.step=  vector of maximal parameter changes per iteration
%
% output:  OPTIMUM is the parameter set that maximizes the instrument output
%            or the MONITORS for single simulations.
%          MONITORS contains the instrument output (e.g. at OPTIMUM).
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% example: 
%   [p,f]=mcstas('templateDIFF',struct('RV',1),struct('mpi',2,'TolFun','0.1%'))
%   subplot(f);
%   monitors=mcstas('templateDIFF',struct('RV',1))
%
% Version: $Revision: 1.2 $
% See also: fminsearch, optimset, http://www.mcstas.org

  if nargin < 2
    error([ 'syntax is: ' mfilename '(instrument, parameters, {options}, {constraints}...)' ] );
  end

  if ~exist('iData')
    error([ mfilename ' requires iFit/iData. Get it at <ifit.mccode.org>. Install it with addpath(genpath(''location/to/iFit''))' ] );
  end
  
  % parse parameter values for mcstas_criteria
  parameter_names = fieldnames(parameters);
  fixed_names     = {};
  variable_names  = {};
  fixed_pars      = {};
  variable_pars   = [];
  for index=1:length(parameter_names)
    value       = getfield(parameters, parameter_names{index});
    if isnumeric(value)
      variable_pars(end+1) = mean(value(:));
      variable_names{end+1}= parameter_names{index};
    else
      fixed_pars{end+1} = char(value);
      fixed_names{end+1}= parameter_names{index};
    end
  end
  
  % McStas configuration and end-user choices
  options.variable_names = variable_names;
  options.fixed_names    = fixed_names;
  options.fixed_pars     = fixed_pars;
  pars                   = variable_pars;
  options.instrument     = instrument;
  if ~isfield(options,'dir'),       options.dir       = tempname;   use_temp=1; else use_temp=0; end
  if ~isfield(options,'ncount'),    options.ncount    = 1e5;        end
  if  isfield(options,'compile') && (options.compile || strcmp(options.compile,'yes'))
    mcstas_criteria(pars, options);
    options = rmfield(options,'compile');
  end
  
  % define simulation or optimization mode (if not set before)
  if ~isfield(options,'mode')
    if isfield(options, 'optimizer') | isfield(options,'TolFun') | isfield(options,'TolX') | isfield(options,'type')
      options.mode      = 'optimize'; 
    else
      options.mode      = 'simulate';
    end
  end
  
  if strcmp(options.mode,'optimize')
    if ~isfield(options,'type'),      options.type      = 'maximize'; end
    if ~isfield(options,'optimizer'), options.optimizer = @fminimfil; end
    
    if isfield(options,'monitors')
      options.monitors = cellstr(options.monitors);
    end
    
    % specific optimizer configuration
    optimizer_options = feval(options.optimizer,'defaults');
    optimizer_options.OutputFcn = 'fminplot';
    optimizer_options.Display   = 'iter';
    field_names=fieldnames(optimizer_options);
    for index=1:length(fieldnames(optimizer_options))
      if ~isfield(options, field_names{index})
        options = setfield(options, field_names{index}, getfield(optimizer_options, field_names{index}));
      end
    end
  
  
    if nargin < 4
    % launch optimizer
    [pars,fval,exitflag,output] = feval(options.optimizer, ...
      @(pars) mcstas_criteria(pars, options), pars, options);
    else
    [pars,fval,exitflag,output] = feval(options.optimizer, ...
      @(pars) mcstas_criteria(pars, options), pars, options, varargin{:});
    end
    options.mode = 'simulate';
    fval = mcstas_criteria(pars, options);
  else
    % single simulation
    fval     = mcstas_criteria(pars, options);
    pars     = fval;
    exitflag = 0;
    output   = options;
  end
  
  if use_temp==1
    % clean up last iteration result when stored into a temporary location
    success = rmdir(options.dir,'s');
  end
  
% end of mcstas function

% ------------------------------------------------------------------------------
function criteria = mcstas_criteria(pars, options)

  % clean up previous iteration result
  success = rmdir(options.dir,'s');
  
  % launch simulation
  cmd = [ 'mcrun ' options.instrument ' --dir=' options.dir ' --ncount=' num2str(options.ncount) ];
  if isfield(options,'gravitation') && (options.gravitation || strcmp(options.gravitation,'yes'))
    cmd = [ cmd ' --gravitation' ];
  end
  if isfield(options,'compile') && (options.compile || strcmp(options.compile,'yes'))
    cmd = [ cmd ' --force-compile' ];
  end
  if isfield(options,'mpi')
    cmd = [ cmd ' --mpi=' num2str(options.mpi) ];
  end
  if isfield(options,'seed')
    cmd = [ cmd ' --seed=' num2str(options.seed) ];
  end
  for index=1:length(options.variable_names)
    cmd = [ cmd ' ' options.variable_names{index} '=' num2str(pars(index)) ];
  end
  for index=1:length(options.fixed_names)
    cmd = [ cmd ' ' options.fixed_names{index} '=' options.fixed_pars{index} ];
  end
  disp(cmd);
  system(cmd);
  
  if nargout ==0, return; end
  directory = options.dir;

  if ~isempty(dir([ directory filesep 'mcstas.sim' ]))
    directory = [ directory filesep 'mcstas.sim' ];
  end

  % import McStas simulation result
  sim = iData(directory); % a vector of monitors (iData objects)
  
  if strcmp(options.mode,'optimize')
    if isfield(options,'monitors') & numel(sim) > 1
      % restrict monitors from simulation by matching patterns
      use_monitors = zeros(size(sim));
      for index=1:length(options.monitors)
        % find monitors that match a search token
        this = cellfun('isempty', findstr(sim, options.monitors{index}));
        this = find(this == 0); % find those that are not empty
        use_monitors(this) = 1;
      end
      if any(use_monitors)
        sim = sim(find(use_monitors));
      end
    end
    criteria = [];
    for index=1:length(sim)
      this = double(sim(index));
      this = sum(this(:));
      if strcmp(options.type,'maximize')
        this = -sum(this);
      end % else minimize
      criteria(end+1) = this;
    end
  else
    criteria = sim;
  end

% end of mcstas_criteria (inline mcstas)
