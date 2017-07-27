function [istop, message] = inline_private_check(pars, fval, funccount, options, constraints)
% standard checks
% inline_private_check(pars, fval, funccount, options, constraints)
% or
% options=inline_private_check(options, default_options);

  istop=0; message='';
  
  % check of option members
  if nargin==2 % check(pars, defaults)
    options=pars;
    if ~isfield(options,'MinFunEvals'), options.MinFunEvals=0; end
    default=fval; % initial values from optimget, as a structure
    checks =fieldnames(default);
    default_static.TolFun     =1e-3;  % these are necessary in any case
    default_static.TolX       =1e-8;
    default_static.MaxIter    =1000;
    default_static.MaxFunEvals=10000;
    default_static.Display    ='';
    default_static.FunValCheck='off';
    default_static.OutputFcn  =[];
    default_static.MinFunEvals=0;
    default_static.PlotFcns   =[];
    % replace required options which are empty
    for f=fieldnames(default_static)'
      if (isfield(options, f{1}) && isempty(options.(f{1}))) || ~isfield(options, f{1})
        options.(f{1}) = default_static.(f{1});
      end
    end
    % put missing options obtained from the optimizer itself
    for f=fieldnames(default)'
      if ~isfield(options, f{1}) options.(f{1}) = default.(f{1}); end
    end
    
    istop=options;
    return
  end
  
  % reduce pars and constraints to only those that vary
  for f={'min','max','step','set'}
      if isfield(constraints, f{1})
          c=constraints.(f{1}); c=c(constraints.index_variable); constraints.(f{1})=c;
      end
  end
  pars = pars(constraints.index_variable);
  
  pars_prev = constraints.parsPrevious;
  fval_prev = constraints.criteriaPrevious;
  fval_best = constraints.criteriaBest;
  fval_mean = mean(constraints.criteriaHistory(max(1, length(constraints.criteriaHistory)-10):end));
  
  % handle relative stop conditions
  if isfield(options,'TolFunChar')
    options.TolFun = options.TolFunChar;
  end
  if ischar(options.TolFun)
    if options.TolFun(end)=='%'
      options.TolFun(end)='';
      options.TolFun = abs(str2num(options.TolFun)*fval/100);
    else
      options.TolFun = str2num(options.TolFun);
    end
  end
  if isfield(options,'TolXChar')
    options.TolX = options.TolXChar;
  end
  if ischar(options.TolX)
    if options.TolX(end)=='%'
      options.TolX(end)='';
      options.TolX = abs(str2num(options.TolX)*pars/100);
    else
      options.TolX = str2num(options.TolX);
    end
  end

  % normal terminations: function tolerance reached
  if isempty(options.MinFunEvals) || funccount >= options.MinFunEvals
    if ~isempty(options.TolFun) && options.TolFun ~= 0 && funccount >= 5*length(pars)
      if (all(0 < fval) && all(fval <= options.TolFun)) % stop on lower threshold
        istop=-1;
        message = [ 'Converged: Termination function tolerance criteria reached (fval <= options.TolFun=' ...
                  num2str(options.TolFun) ')' ];
      end
      if ~istop
        % stop on criteria change
        if  all(abs(fval-fval_prev) < options.TolFun) ...
         && all(abs(fval-fval_prev) > 0) ...
         && all(fval < fval_mean - options.TolFun) 
          istop=-12;
          message = [ 'Converged: Termination function change tolerance criteria reached (delta(fval) < options.TolFun=' ...
                  num2str(options.TolFun) ')' ];
        end
      end
    end
  end
  
  % normal terminations: parameter variation tolerance reached, when function termination is also true
  if (istop==-1 || istop==-12) 
    if ~isempty(options.TolX) && isnumeric(options.TolX)
      index=find(isfinite(options.TolX) & options.TolX);
      if all(abs(pars(index)-pars_prev(index)) < abs(options.TolX(index))) ...
      && any(abs(pars(index)-pars_prev(index)) > 0)
        istop=-5;
        message = [ 'Converged: Termination parameter tolerance criteria reached (delta(parameters) <= options.TolX=' ...
              num2str(mean(options.TolX)) ')' ];
      end
    end
  end

  % abnormal terminations
  if options.MaxFunEvals > 0 & funccount >= options.MaxFunEvals
    istop=-3;
    message = [ 'Maximum number of function evaluations reached (options.MaxFunEvals=' ...
              num2str(options.MaxFunEvals) ')' ];
  end

  % the function value is nan or parameters just went to nan
  if strcmp(options.FunValCheck,'on') && (any(isnan(fval) | isinf(fval)))
    istop=-4;
    message = 'Function value is Inf or Nan (options.FunValCheck)';
  end
  
  % make sure parameters remain non NaN's for e.g. plotting
  if any(isnan(pars))
    index = find(isnan(pars(:)) & ~isnan(pars_prev(:)));
    if ~isempty(index)
      pars(index) = pars_prev(index);
    end
    index = find(isnan(pars(:)) & ~isnan(constraints.parsBest(:)));
    if ~isempty(index)
      pars(index) = constraints.parsBest(index);
    end
    index = find(isnan(pars(:)) & ~isnan(constraints.parsStart(:)));
    if ~isempty(index)
      pars(index) = constraints.parsStart(index);
    end
  end
    
  % create an optimValues structure compliant with OutputFcn and PlotFcns
  optimValues = options;
  if ~isfield(optimValues,'state')
    if istop,               optimValues.state='done';
    else                    optimValues.state='iter'; end
  end
  optimValues.funcount   = funccount;
  optimValues.funcCount  = funccount;
  optimValues.funccount  = funccount;
  optimValues.iteration  = funccount;
  optimValues.fval       = sum(fval(:));
  if isfield(options,'procedure'),        optimValues.procedure=options.procedure;
  elseif isfield(options, 'algorithm'),   optimValues.procedure=options.algorithm;
  elseif isfield(options, 'optimizer'),   optimValues.procedure=options.optimizer;
  else optimValues.procedure  = 'iteration'; end
    
  % we assemble the list of external functions to call: PlotFcns OutputFcn
  % syntax: stop = OutputFcn(pars, optimValues, state)
  if ~isempty(options.PlotFcns) && ~iscell(options.PlotFcns)
    ExternalFcns = { options.PlotFcns }; 
  else ExternalFcns = options.PlotFcns; end
  if ~isempty(options.OutputFcn), ExternalFcns{end+1} = options.OutputFcn; end

  for index=1:numel(ExternalFcns)
    PlotFcns = ExternalFcns{index};
    if ischar(PlotFcns) || isa(PlotFcns, 'function_handle')
      try
        istop2 = feval(PlotFcns, pars, optimValues, optimValues.state);
      catch ME
        disp(getReport(ME))
        istop2 = 0; % failed ExternalFcns ignored
      end
      if istop2 && ~istop
        istop=-6;
        message = [ 'Algorithm was terminated by the output function ' ...
                    inline_localChar(PlotFcns) ];
        break
      end
    end
  end % ExternalFcns
    
  if istop
    funccount = -funccount; % trigger iteration display
  end

end % inline_private_check
