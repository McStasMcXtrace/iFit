function [pars,fval,exitflag,output] = mcstas(instrument, parameters, options)
% [OPTIMUM,MONITORS,EXITFLAG,OUTPUT] = mcstas(INSTRUMENT, PARAMETERS, OPTIONS) : run and optimize a McStas simulation
%
% A wrapper to the McStas package to either execute a simulation, or optimize a
%   set of parameters. To select the optimization mode, set the options.mode='optimize'
%   or any other optimization configuration parameter (TolFun, TolX, ...). Default 
%   mode is 'simulate', which also includes scanning capability.
% The syntax:
%   mcstas(instrument,'compile')     assemble and compile the instrument
%   mcstas(instrument,'compile mpi') same with MPI support
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
%           parameters.name = string (simulation and optimization)
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
%             the monitor names can contain expressions made of the monitor name 
%             followed by any expression using 'this' to refer to the monitor content
%             such as in 'Monitor1/std(this,1)' which divides the Signal by its X peak width.
%           options.help:   set it to 'yes' or 1 to get help on the instrument and exit
%           options.info:   set it to 'yes' or 1 to get information on the instrument and exit
%           options.optimizer: function name of the optimizer to use (string or function handle)
%           options.OutputFcn: monitors the scan/optimization process on a plot (string)
%         as well as other optimizer options such as
%           options.TolFun ='0.1%'   stop when criteria changes are smaller that 0.1%
%           options.Display='final'
%
% output:  OPTIMUM is the parameter set that maximizes the instrument output, or
%            the integral monitor values for the simulation (as iData object)
%          MONITORS contains the instrument output as iData objects. Each object has an additional
%            Parameter member alias which holds the instrument parameters.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% example: to optimize RV, display result, and then perform a scan
%   [p,f]=mcstas('templateDIFF', struct('RV',[0.5 1 1.5]), struct('TolFun','0.1%','monitors','Banana'))
%   subplot(f); disp(f.Parameters)
%   [monitors_integral,scan]=mcstas('templateDIFF' ,struct('RV',[0.5 1 1.5]))
%   plot(monitors_integral)
%
% Version: $Revision: 1.22 $
% See also: fminsearch, fminimfil, optimset, http://www.mcstas.org

% inline: mcstas_criteria

  if nargin < 1
    error([ 'syntax is: ' mfilename '(instrument, parameters, {options})' ] );
  end

  if ~exist('iData')
    error([ mfilename ' requires iFit/iData. Get it at <ifit.mccode.org>. Install it with addpath(genpath(''/path/to/iFit''))' ] );
  end
  
% PARSE INPUT ARGUMENTS ========================================================
  
  if nargin > 2 && ischar(options)
    options= str2struct(options);
  end
  
  options.instrument     = instrument;
  % define simulation or optimization mode (if not set before)
  if ~isfield(options,'mode')
    if isfield(options, 'optimizer') || isfield(options,'TolFun') | ...
       isfield(options,'TolX') | isfield(options,'type') | ...
       isfield(options,'MaxFunEvals') | isfield(options,'MaxIter')
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
      if ~isfield(options,'TolFun')
        options.TolFun = '0.1%';
      end
    else
      options.ncount    = 1e6;
    end
  end
  % force compile before going further ?
  if  isfield(options,'compile') & (options.compile | strcmp(options.compile,'yes'))
    ncount = options.ncount;
    options.ncount=0;
    mcstas_criteria([], options);
    options = rmfield(options,'compile');
    options.ncount = ncount;
  end  
  
  % parse parameter values for mcstas_criteria
  if ischar(parameters) && ~isempty(strfind(parameters,'compile'))
    options.compile = 1;
    if ischar(parameters) && (~isempty(strfind(parameters,' mpi')) || ~isempty(strfind(parameters,'mpi ')))
      options.mpi=2;
    end
    parameters      = [];
  end
  if nargin < 1, parameters = []; end
  if ~isempty(parameters) && ischar(parameters)
    parameters = str2struct(parameters);
  end
  if ~isempty(parameters) && isstruct(parameters)
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
  end % for index
  
  % optimizer configuration and end-user choices
  options.variable_names = variable_names;
  options.variable_pars  = variable_pars;
  options.fixed_names    = fixed_names;
  options.fixed_pars     = fixed_pars;
  options.scan_size      = scan_size;
  pars                   = variable_pars;
  
% Launch optimize and simulate mode ============================================

  if isfield(options,'monitors')
    options.monitors = cellstr(options.monitors);
  end

  if strcmp(options.mode,'optimize') % ================================ OPTIMIZE
    % optimize simulation parameters
    if ~isfield(options,'type'),      options.type      = 'maximize'; end
    if ~isfield(options,'optimizer'), options.optimizer = @fminpso; end
    
    % specific optimizer configuration
    optimizer_options = feval(options.optimizer,'defaults');
    optimizer_options.OutputFcn = 'fminplot';
    optimizer_options.TolFun    = '1%';
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
    output.parameters=parameters;
    options.mode = 'simulate'; % evaluate best solution, when optimization is over
    if nargout > 1
      fprintf(1,'Evaluating final solution...\n');
      [dummy,fval] = mcstas_criteria(pars, options);
      output.command=get(fval(1), 'Execute');
    end
    % re-create a structure of parameters
    pars_struct = [];
    for index=1:length(pars)
      pars_struct.(variable_names{index}) = pars(index);
    end
    pars = pars_struct;
  else % ================================================ single simulation/scan

    [p, fval] = mcstas_criteria(pars, options);   % may fail at execution
    
    fval = squeeze(fval);
    p    = squeeze(p);
    if iscell(fval) && ~isempty(fval)
      % before converting to a single iData array, we check that all
      % simulations returned the same number of monitors
      for index=1:numel(fval)
          if isempty(fval{index}), fval{index}=iData; end
      end
      siz = cellfun('prodofsize',fval);
      if all(siz == siz(1))
        fval=iData(fval);
      end
    end
    if isempty(fval)
      pars=[]; fval=[];
    elseif nargout < 2
      pars     = fval;
    else
      a = iData(p);
      try
        t = fval(1); isscan = t.Data.Scan;
        t = ' Scan of'; isscan = 1;
      catch
        t = ''; isscan = 0;
      end
      for index=1:length(options.variable_names)
        setalias(a, options.variable_names{index}, options.variable_pars{index});
        if isscan==1
          setaxis(a, index, options.variable_names{index});
        end
        t = [ t ' ' options.variable_names{index} ];
      end
      setalias(a, 'Criteria', 1:size(a, ndims(a)), 'Monitor index');
      if isscan==1
        setaxis(a, length(options.variable_names)+1, 'Criteria');
      end
      % add other metadata
      set(a, 'Data.Parameters', parameters);
      set(a, 'Data.Criteria', p);
      if iscell(fval)
      set(a, 'Data.Execute', get(fval{1},'Execute'));
      else
      set(a, 'Data.Execute', get(fval(1),'Execute'));
      end
      set(a, 'Data.Options', options);
      setalias(a, 'Parameters', 'Data.Parameters','Instrument parameters');
      setalias(a, 'Execute', 'Data.Execute','Command line used for Mcstas execution');
      setalias(a, 'Options', 'Data.Options','Options used for Mcstas execution');
      a.Title = [ instrument ':' t ];
      a.Label = instrument;
      pars = a;
    end 
    exitflag = 0;
    output   = options;
  end % else single simulation
  
  if use_temp==1
    % clean up last iteration result when stored into a temporary location
    success = rmdir(options.dir,'s');
  end

end
% end of mcstas function

% ------------------------------------------------------------------------------

function system_wait(cmd, options)
% inline function to execute command and wait for its completion (under Windows)
% dots are displayed under Windows after every minute waiting.
  [status, result]=system(cmd);
  disp(result);
  if status ~= 0, return; end
  if ispc % wait for completion by monitoring the number of elements in the result directory
    t=tic; t0=t; first=1;
    a=dir(options.dir);
    while length(a) <= 3 % only 'mcstas.sim', '.', '..' when simulation is not completed yet
      if toc(t) > 60
        if first==1 % display initial waiting message when computation lasts more than a minute
            fprintf(1, 'mcstas: Waiting for completion of %s simulation (dots=minutes).\n', options.instrument);
        end
        fprintf(1,'.');
        t=tic; first=first+1;
        if first>74 % go to next line when more than 75 dots in a row...
            first=2;
            fprintf(1,'\n');
        end
      end
      a=dir(options.dir); % update directory content list
    end
    fprintf(1,' DONE [%10.2g min]\n', toc(t0)/60);
  end
  % wait for directory to be 'stable' (not being written)
  this_sum=0;
  while 1
      a=dir(options.dir);
      new_sum=0;
      for index=1:length(a)
          new_sum = new_sum+a(index).datenum;
      end
      if this_sum == new_sum, break; end
      this_sum = new_sum;
      pause(1);
  end
end % system_wait

% ------------------------------------------------------------------------------

function [criteria, sim, ind] = mcstas_criteria(pars, options, criteria, sim, ind)
% inline function to compute a single simulation, or a vector of simulations (recursive calls)
  
  % launch simulation with mcrun
  cmd = [ 'mcrun ' options.instrument ];
  
  % usual McStas/mcrun options
  if isfield(options,'compile') & (options.compile | strcmp(options.compile,'yes'))
    cmd = [ cmd ' --force-compile' ];
    if isempty(pars)
      options.ncount=0;
      options.dir   ='';
    end
  end
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
  if isfield(options,'mpi')
    if isempty(options.mpi)
      cmd = [ cmd ' --mpi' ];
    elseif options.mpi > 1
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
      if isnumeric(pars) % all numerics
        cmd = [ cmd ' ' options.variable_names{index} '=' num2str(pars(index)) ];
      else % some are cells and chars
        % scan mode, with vector parameters
        this = pars{index};
        if isnumeric(this) && length(this) == 1
          cmd = [ cmd ' ' options.variable_names{index} '=' num2str(this) ];
        elseif ischar(this)
          cmd = [ cmd ' ' options.variable_names{index} '=' this ];
        elseif isvector(this) % parameter is a vector of numerics/scans
          for index_pars=1:length(this) % scan the vector parameter elements
            ind{index} = index_pars; % coordinates of this scan step in the parameter space indices
            if isnumeric(this)
              pars{index} = this(index_pars);
            elseif iscell(this)
              pars{index} = this{index_pars};
            end
            % recursive call to handle all scanned parameters
            [this_criteria, this_sim, ind] = mcstas_criteria(pars, options, criteria, sim, ind);
            if ~iscell(this_sim)   % single simulation
              if isempty(criteria) % initialize arrays to the right dimension
                criteria = zeros([ options.scan_size length(this_criteria) ]); % array of 0
                sim      = cell( size(criteria) );                             % empty cell
              end

              % add single simulation to scan arrays
              % store into the last dimensionality (which holds monitors and integrated values)
              if ~isempty(this_sim)
              for index_mon=1:length(this_criteria)
                if isempty(ind), this_ind = { index_mon };
                else this_ind = { ind{:} index_mon }; end
                this_ind(cellfun('isempty',this_ind))={1};
                try
                sim{      sub2ind(size(sim), this_ind{:}) } = this_sim(index_mon);
                criteria( sub2ind(size(sim), this_ind{:}) ) = this_criteria(index_mon);
                catch
                sim{      this_ind{:} } = this_sim(index_mon);
                criteria( this_ind{:})  = this_criteria(index_mon);
                end
                
              end
              end
            else
              criteria = this_criteria;
              sim      = this_sim;
            end
            % optionally plot the criteria during the scan...
              if isfield(options, 'OutputFcn') 
                if ~isempty(options.OutputFcn)
                  if isvector(criteria)
                    plot(criteria)
                    xlabel('Scan step'); ylabel('Monitors'); 
                    t=title([ options.instrument ': ' options.variable_names{index} '=' num2str(pars{index}) ]); 
                    set(t,'interpreter','none');
                    drawnow
                  elseif length(size(criteria)) == 2
                    surf(criteria);
                    t=title([ options.instrument ': ' options.variable_names{index} '=' num2str(pars{index}) ]); 
                    set(t,'interpreter','none');
                    drawnow
                  end
                end
              end
          end % for index_pars
          if nargout < 2
            criteria = sim;
          end
          return  % return from scan
        end % elseif isvector(this): parameter value given as vector 
      end
    end % for index
  end
  % non scanned parameters (chars, fixed)
  if isfield(options,'fixed_names')
    for index=1:length(options.fixed_names)
      cmd = [ cmd ' ' options.fixed_names{index} '=' options.fixed_pars{index} ];
    end
  end
  
  % Execute simulation =========================================================
  disp([ mfilename ': ' options.mode ' ' cmd ]);
  system_wait(cmd, options);
  
  if nargout ==0, return; end
  if isfield(options,'ncount') && options.ncount == 0
    return
  end
  directory = options.dir;

  % import McStas simulation result
  sim = [];
  try
    % first try to import monitors from their file names
    if isfield(options,'monitors')
      for index=1:length(options.monitors)
        [name, R] = strtok(options.monitors{index},' ,;/*+-(){}:%$.');
        sim = [ sim iData(fullfile(directory,[ name '*' ])) ];
        setalias(sim, 'CriteriaExpression', R);
      end
    end
    if isempty(sim)
      % if designated monitor file name import fails, import all simulation content
      if ~isempty(dir([ directory filesep 'mcstas.sim' ]))
        directory = [ directory filesep 'mcstas.sim' ];
      end
      sim = iData(directory); % a vector of monitors (iData objects)
      
      % filter all simulation monitor
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
    end
  catch
    criteria=0; sim=[]; ind=[];
    return
  end

  
  % option to plot the monitors
  if isfield(options, 'OutputFcn')
    % is this window already opened ?
    h = findall(0, 'Tag', 'mcstasplot');
    if isempty(h) % create it
      h = figure('Tag','mcstasplot', 'Unit','pixels');
      tmp = get(h, 'Position'); tmp(3:4) = [500 400];
      set(h, 'Position', tmp);
    end

    % raise existing figure (or keep it hidden) and add parameters on top
    if gcf ~= h, figure(h); end
    hold off
    subplot(sim,'view2 axis tight');
    ud.Parameters = get(sim(1),'Parameters');
    ud.Execute=cmd;
    ud.pars   =pars;
    ud.Options=options;
    set(h, 'UserData', ud);
    xl=xlim; yl=ylim;
    f = fieldnames(ud.Parameters);
    c = struct2cell(ud.Parameters);
    if length(f) > 20, f=f(1:20); c=c(1:20); end
    s = class2str('p',cell2struct(c(:),f(:),1),'no comment');
    text(xl(1), mean(yl), s,'Interpreter','none');
    hold off
  end
  
  % evaluate the criteria specifications (monitors, expressions)
  criteria = zeros(length(sim),1);
  for index=1:length(sim)
    R = '';
    try
      R = getalias(sim(index), 'CriteriaExpression');
      this = sim(index);
      eval([ 'this = this' R ';' ]);
    end
    this = double(this);
    this = sum(this(:));
    if isfield(options,'type') & strcmp(options.type,'maximize')
      this = -sum(this);
    end % else minimize
    criteria(index) = this;
  end
  
  % add aliases to the output objects (Parameters, Command line, ...)
  if nargout > 1
    if isnumeric(pars)
      this_pars = cell(size(pars));
      for index=1:length(this_pars)
        this_pars{index} = pars(index);
      end
    else
      this_pars = pars;
    end
    c = { this_pars{:} , options.fixed_pars{:} }; c=c(:);
    f = { options.variable_names{:} , options.fixed_names{:}}; f=f(:);
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
end % mcstas_criteria
% end of mcstas_criteria (inline mcstas)
