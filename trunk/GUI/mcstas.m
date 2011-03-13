function [pars,fval,exitflag,output] = mcstas(instrument, parameters, options)
% [OPTIMUM,MONITORS,EXITFLAG,OUTPUT] = mcstas(INSTRUMENT, PARAMETERS, OPTIONS) : run and optimize a McStas simulation
%
% A wrapper to the McStas package to either execute a simulation, or optimize a
%   set of parameters. To select the optimization mode, set the options.mode='optimize'
%   or any other optimization configuration parameter (TolFun, TolX, ...). Default 
%   mode is 'simulate', which also includes scanning capability.
%
% input:  INSTRUMENT: name of the instrument description to run (string)
%         PARAMETERS: a structure that gives instrument parameter names and values (structure)
%           parameters.name = scalar (optimization and simulation)
%             define a single value      when options.mode='simulation'
%             define the starting value  when options.mode='optimization'
%           parameters.name = vector and cell array (simulation)
%             define a scanning range for a series of simulations, when options.mode='simulation'
%             multi-dimensional scans (that is more than one parameter given as 
%             vector) are possible.
%           parameters.name = vector (optimize)
%             the vector can be 3 elements to define [min start max]. 
%             the vector can be 2 elements to define [min max]. 
%             the vector can be 4+ elements to define [min ... max]. 
%             For 2 and 4+ vectors, the parameter starting value is the mean of the vector.
%           parameters.name = string (simulation ans optimization)
%             defines a fixed parameter, that can not be optimized.
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
%           options.help:   set it to 'yes' or 1 to get help on the instrument and exit
%           options.info:   set it to 'yes' or 1 to get information on the instrument and exit
%           options.optimizer: function name of the optimizer to use (string or function handle)
%         as well as other optimizer options such as
%           options.TolFun ='0.1%'   stop when criteria changes are smaller that 0.1%
%           options.Display='final'
%
% output:  OPTIMUM is the parameter set that maximizes the instrument output, or
%            the integral monitor values for the simulation (as idata object)
%          MONITORS contains the instrument output as iData objects.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% example: to optimize RV, display result, and then perform a scan
%   [p,f]=mcstas('templateDIFF', struct('RV',1), struct('TolFun','0.1%'))
%   subplot(f);
%   [monitors_integral,scan]=mcstas('templateDIFF' ,struct('RV',[0.5 1 1.5]))
%   plot(monitors_integral)
%
% Version: $Revision: 1.4 $
% See also: fminsearch, fminimfil, optimset, http://www.mcstas.org

% inline: mcstas_criteria

  if nargin < 1
    error([ 'syntax is: ' mfilename '(instrument, parameters, {options})' ] );
  end

  if ~exist('iData')
    error([ mfilename ' requires iFit/iData. Get it at <ifit.mccode.org>. Install it with addpath(genpath(''location/to/iFit''))' ] );
  end
  
% PARSE INPUT PARAMETERS =======================================================
  
  options.instrument     = instrument;
  % define simulation or optimization mode (if not set before)
  if ~isfield(options,'mode')
    if isfield(options, 'optimizer') | isfield(options,'TolFun') | ...
       isfield(options,'TolX') | isfield(options,'type') 
      options.mode      = 'optimize'; 
    else
      options.mode      = 'simulate';
    end
  end
  
  if ~isfield(options,'dir')
    options.dir = tempname;
    use_temp    = 1; 
  else 
    use_temp=0; 
  end
  if ~isfield(options,'ncount')
    if strcmp(options.mode, 'optimize')
      options.ncount    = 1e5;
    else
      options.ncount    = 1e6;
    end
  end
  % force compile before going further ?
  if  isfield(options,'compile') & (options.compile | strcmp(options.compile,'yes'))
    ncount = options.ncount;
    options.ncount=0
    mcstas_criteria([], options);
    options = rmfield(options,'compile');
    options.ncount = ncount;
  end  
  
  % parse parameter values for mcstas_criteria
  if nargin < 1, parameters = []; end
  if ~isempty(parameters)
    parameter_names = fieldnames(parameters);
  else
    parameter_names = {};
  end  
  fixed_names     = {};
  variable_names  = {};
  fixed_pars      = {};
  scan_size       = [];
  if strcmp(options.mode, 'optimize')
    variable_pars   = [];
  else
    variable_pars   = {};
  end
  constraints = []; % for optimization min/max
  for index=1:length(parameter_names)
    value = getfield(parameters, parameter_names{index});
    value = value(:);
    if ischar(value) % fixed values for both simulation and optimization
      fixed_pars{end+1} = value(:)';
      fixed_names{end+1}= parameter_names{index};
    elseif isvector(value) % this is a vector: numeric or cell
      if strcmp(options.mode, 'optimize')
        if ~isnumeric(value)
          error([ mfilename ': Parameter ' parameter_names{index} ' of type ' class(value) ' is not supported.' ...
            sprintf('\n') 'Optimization only supports numerics as input.']);
        end
        % optimize parameter=numerical vector 
        if length(value) > 1
          % the vector in optimization mode allows definition of constraints
          if isempty(constraints) % define constraints to NaN up to now
            constraints.min = NaN*ones(size(variable_pars));
            constraints.max = constraints.min;
          end
          constraints.min(end+1) = min(value); % add the new constraints for this parameter.
          constraints.max(end+1) = max(value);
          if length(value) == 3
            value = value(2); % syntax: parameter=[ min start max ]
          end
          variable_pars(end+1) = mean(value(:));
        else
          if ~isempty(constraints)
            constraints.min(end+1) = NaN; % if use constraints elsewhere,
            constraints.max(end+1) = NaN; % fill in no-constraints for this parameter
          end
          variable_pars(end+1) = value;
        end
      else % simulate mode: we just store the vector, that will be scanned through
        scan_size(end+1)     = length(value);
        variable_pars{end+1} = value;
      end
      variable_names{end+1}= parameter_names{index};
    else % not a vector (probably a char)
      error([ mfilename ': Parameter ' parameter_names{index} ' of type ' class(value) ' is not supported.' ...
        sprintf('\n') 'Prefer numerical or cell vectors.']);
    end
  end
  
  % optimizer configuration and end-user choices
  options.variable_names = variable_names;
  options.variable_pars  = variable_pars;
  options.fixed_names    = fixed_names;
  options.fixed_pars     = fixed_pars;
  options.scan_size      = scan_size;
  pars                   = variable_pars;
  
% Launch optimize and simulate mode ============================================

  if strcmp(options.mode,'optimize')
    if ~isfield(options,'type'),      options.type      = 'maximize'; end
    if ~isfield(options,'optimizer'), options.optimizer = @fminimfil; end
    
    if isfield(options,'monitors')
      options.monitors = cellstr(options.monitors);
    end
    
    % specific optimizer configuration
    optimizer_options = feval(options.optimizer,'defaults');
    optimizer_options.OutputFcn = 'fminplot';
    optimizer_options.TolFun    = '0.1%';
    field_names=fieldnames(optimizer_options);
    for index=1:length(fieldnames(optimizer_options))
      if ~isfield(options, field_names{index})
        options = setfield(options, field_names{index}, getfield(optimizer_options, field_names{index}));
      end
    end

    % launch optimizer with or without constraints
    if isempty(constraints)
      [pars,fval,exitflag,output] = feval(options.optimizer, ...
        @(pars) mcstas_criteria(pars, options), pars, options);
    else
      [pars,fval,exitflag,output] = feval(options.optimizer, ...
        @(pars) mcstas_criteria(pars, options), pars, options, constraints);
    end
    options.mode = 'simulate'; % evaluate best solution, when optimization is over
    if nargout > 1
      [dummy,fval] = mcstas_criteria(pars, options);
    end
  else
    % single simulation/scan
    [pars, fval] = mcstas_criteria(pars, options);
    if iscell(fval)
      % before converting to a single iData array, we check that all
      % simulations returned the same number of monitors
      siz = cellfun('prodofsize',fval);
      if all(siz == siz(1))
        siz = [ siz(1) size(fval) ]; % size of the array to generate
        fval = reshape(iData(fval), siz);
        % permute dimensions in order to have monitors as last
        perm = 1:length(siz); perm(1) = perm(end); perm(end) = 1;
        fval = permute(fval, perm);
      end
    end
    if nargout < 2
      pars     = fval;
    else
      a = iData(pars);
      t = 'Scan of';
      for index=1:length(options.variable_names)
        setalias(a, options.variable_names{index}, options.variable_pars{index});
        setaxis(a, index, options.variable_names{index});
        t = [ t ' ' options.variable_names{index} ];
      end
      setalias(a, 'Criteria', 1:size(a, ndims(a)), 'Monitor index');
      setaxis(a, length(options.variable_names)+1, 'Criteria');
      a.Title = [ instrument ': ' t ];
      a.Label = instrument;
      pars = a;
    end
    exitflag = 0;
    output   = options;
  end
  
  if use_temp==1
    % clean up last iteration result when stored into a temporary location
    success = rmdir(options.dir,'s');
  end

end
% end of mcstas function

% ------------------------------------------------------------------------------

function [criteria, sim, ind] = mcstas_criteria(pars, options, criteria, sim, ind)
% inline function to compute a single simulation, or a vector of simulations (recursive calls)
  
  % launch simulation with mcrun
  cmd = [ 'mcrun ' options.instrument ];
  
  % usual McStas/mcrun options
  if isfield(options,'ncount') && ~isempty(options.ncount)
    cmd = [ cmd ' --ncount=' num2str(options.ncount) ];
  end
  if isfield(options,'dir') && ~isempty(options.dir)
    % clean up previous simulation result
    success = rmdir(options.dir,'s');
    cmd = [ cmd ' --dir=' options.dir ];
  end
  if isfield(options,'gravitation') && (options.gravitation || strcmp(options.gravitation,'yes'))
    cmd = [ cmd ' --gravitation' ];
  end
  if isfield(options,'compile') & (options.compile | strcmp(options.compile,'yes'))
    cmd = [ cmd ' --force-compile' ];
  end
  if isfield(options,'mpi')
    if isempty(options.mpi)
      cmd = [ cmd ' --mpi' ];
    else
      cmd = [ cmd ' --mpi=' num2str(options.mpi) ];
    end
  end
  if isfield(options,'info')
    cmd = [ cmd ' --info' ];
  end
  if isfield(options,'help')
    cmd = [ cmd ' --help' ];
  end
  if isfield(options,'seed') && ~isempty(options.seed)
    cmd = [ cmd ' --seed=' num2str(options.seed) ];
  end
  
  % handle single simulation and vectorial scans
  if isfield(options,'variable_names')
    if nargin < 3, 
      ind = cell(1,length(options.variable_names)); ind{1}=1; 
      if ~exist('criteria')
        criteria = [];
        sim      = {};
      end
    end
    for index=1:length(options.variable_names)
      if isnumeric(pars)
        cmd = [ cmd ' ' options.variable_names{index} '=' num2str(pars(index)) ];
      else
        % scan mode, with vector parameters
        this = pars{index};
        if isnumeric(this) && length(this) == 1
          cmd = [ cmd ' ' options.variable_names{index} '=' num2str(this) ];
        elseif ischar(this)
          cmd = [ cmd ' ' options.variable_names{index} '=' this ];
        elseif isvector(this) % parameter is a vector of numerics/scans
          for index_pars=1:length(this) % scan the vector parameter elements
            ind{index} = index_pars;
            if isnumeric(this)
              pars{index} = this(index_pars);
            elseif iscell(this)
              pars{index} = this{index_pars};
            end
            % recursive call to handle all scanned parameters
            [this_criteria, this_sim, ind] = mcstas_criteria(pars, options, criteria, sim, ind);
            if ~iscell(this_sim)
              S.type = '()';
              if isempty(criteria) % initialize arrays to the right dimension
                criteria = zeros([ options.scan_size length(this_sim) ]);
                sim      = cell(options.scan_size);
              end
              % reshape criteria to the last dimensionality
              this_criteria = reshape(this_criteria, [ ones(size(options.scan_size)) length(this_criteria) ]);
              this_sim = reshape(this_sim, [ ones(size(options.scan_size)) length(this_sim) ]);
              % add single simulation to scan arrays
              S.subs   = ind;
              sim      = subsasgn(sim,      S, { this_sim });
              S.subs   = { ind{:}, ':' };
              criteria = subsasgn(criteria, S, this_criteria(:)');
            else
              criteria = this_criteria;
              sim      = this_sim;
            end
          end
          if nargout < 2
            criteria = sim;
          end
          return  % return from scan
        end
      end
    end
  end
  % non scanned parameters (chars, fixed)
  if isfield(options,'fixed_names')
    for index=1:length(options.fixed_names)
      cmd = [ cmd ' ' options.fixed_names{index} '=' options.fixed_pars{index} ];
    end
  end
  disp(cmd);
  system(cmd);
  
  if nargout ==0, return; end
  directory = options.dir;

  if ~isempty(dir([ directory filesep 'mcstas.sim' ]))
    directory = [ directory filesep 'mcstas.sim' ];
  end

  % import McStas simulation result
  try
    sim = iData(directory); % a vector of monitors (iData objects)
  catch
    criteria=0; sim=[]; ind=[];
    return
  end
  
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
  criteria = zeros(length(sim),1);
  for index=1:length(sim)
    this = double(sim(index));
    this = sum(this(:));
    if isfield(options,'type') & strcmp(options.type,'maximize')
      this = -sum(this);
    end % else minimize
    criteria(index) = this;
  end
  if nargout > 1
    if isnumeric(pars)
      this_pars = cell(size(pars));
      for index=1:length(this_pars)
        this_pars{index} = pars(index);
      end
    else
      this_pars = pars;
    end
    c = { this_pars{:} ; options.fixed_pars{:} }; c=c(:);
    f = { options.variable_names{:} ; options.fixed_names{:}}; f=f(:);
    this_pars = cell2struct(c, f, 1);
    set(sim, 'Data.Parameters', this_pars);
    set(sim, 'Data.Criteria', criteria);
    set(sim, 'Data.Execute', cmd);
    set(sim, 'Data.Options', options);
    for index=1:length(f)
      setalias(sim, f{index}, [ 'Data.Parameters.' f{index} ], [ 'Instrument parameter ' f{index} ]);
    end
    setalias(sim, 'Parameters', 'Data.Parameters','Instrument parameters');
    setalias(sim, 'Criteria', 'Data.Criteria','Integral of monitors (criteria)');
    setalias(sim, 'Execute', 'Data.Execute','Command line used for Mcstas execution');
    setalias(sim, 'Options', 'Data.Options','Options used for Mcstas execution');
  end
end
% end of mcstas_criteria (inline mcstas)
