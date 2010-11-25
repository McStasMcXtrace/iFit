function [istop, message] = fmin_private_std_check(pars, fval, iterations, funccount, options, pars_prev, fval_prev)
% standard checks
% fmin_private_std_check(pars, fval, iterations, funccount, options
% or
% fmin_private_std_check(pars, fval, iterations, funccount, options, pars_prev)
% or
% fmin_private_std_check(pars, fval, iterations, funccount, options, pars_prev, fval_prev)
% or
% options=fmin_private_std_check(options);
% or
% options=fmin_private_std_check(options, default_options);


  istop=0; message='';
  
  % check of option members
  if nargin<=2
    options=pars;
    if nargin ==2, 
      default=fval; 
      checks=fieldnames(default);
    else 
      fval=[]; 
      checks={'TolFun','TolX','Display','MaxIter','MaxFunEvals','FunValCheck','OutputFcn','algorithm'};
    end
    
    for index=1:length(checks)
      if ~isfield(options, checks{index}), 
        if isfield(default, checks{index}), 
          options=setfield(options,checks{index},getfield(default, checks{index}));
        else
          options=setfield(options,checks{index},[]); 
        end
      end
    end
    istop=options;
    return
  end

  % normal terminations: function tolerance reached
  if ~isempty(options.TolFun) && options.TolFun
    if (fval <= options.TolFun) & nargin < 7 % stop on lower threshold
      istop=-1;
      message = [ 'Termination function tolerance criteria reached (fval <= options.TolFun=' ...
                num2str(options.TolFun) ')' ];
    end
    if ~istop & nargin >= 7
      if abs(fval-fval_prev) < options.TolFun*abs(fval) ...
       & abs(fval-fval_prev) > 0
        istop=-12;
        message = [ 'Termination function change tolerance criteria reached (delta(fval) < options.TolFun=' ...
                num2str(options.TolFun) ')' ];
      end
    end
  end
  
  % normal terminations: parameter variation tolerance reached, when function termination is also true
  if (istop==-1 || istop==-12) & nargin >= 6
    if ~isempty(options.TolFun) & options.TolX > 0 ...
      & all(abs(pars(:)-pars_prev(:)) < abs(options.TolX*pars(:))) ...
      & any(abs(pars(:)-pars_prev(:)) > 0)
      istop=-5;
      message = [ 'Termination parameter tolerance criteria reached (delta(parameters)/parameters <= options.TolX=' ...
            num2str(options.TolX) ')' ];
    end
  end
  
  % abnormal terminations
  if ~istop
    if options.MaxIter > 0 & iterations >= options.MaxIter
      istop=-2;
      message = [ 'Maximum number of iterations reached (options.MaxIter=' ...
                num2str(options.MaxIter) ')' ];
    end

    if options.MaxFunEvals > 0 & funccount >= options.MaxFunEvals
      istop=-3;
      message = [ 'Maximum number of function evaluations reached (options.MaxFunEvals=' ...
                num2str(options.MaxFunEvals) ')' ];
    end
    
    if strcmp(options.FunValCheck,'on') & any(isnan(fval) | isinf(fval))
      istop=-4;
      message = 'Function value is Inf or Nan (options.FunValCheck)';
    end

    if ~isempty(options.OutputFcn)
      optimValues = options;
      if ~isfield(optimValues,'state')
        if istop,               optimValues.state='done';
        elseif iterations <= 1, optimValues.state='init';
        else                    optimValues.state='iter'; end
      end
      optimValues.iteration  = iterations;
      optimValues.funcount   = funccount;
      optimValues.fval       = fval;
      if isfield(options,'procedure'),        optimValues.procedure=options.procedure;
      elseif isfield(options, 'algorithm'),   optimValues.procedure=options.algorithm;
      else optimValues.procedure  = 'iteration'; end
      istop2 = feval(options.OutputFcn, pars, optimValues, optimValues.state);
      if istop2, 
        istop=-6;
        message = 'Algorithm was terminated by the output function (options.OutputFcn)';
      end
    end
  end
  
% return code     message
%  0                Algorithm terminated normally
% -1                Termination function tolerance criteria reached
% -2                Maximum number of iterations reached
% -3                Maximum number of function evaluations reached
% -4                Function value is Inf or Nan
% -5                Termination parameter tolerance criteria reached
% -6                Algorithm was terminated by the output function
% -7                Maximum consecutive rejections exceeded (anneal)
% -8                Minimum temperature reached (anneal)
% -9                Global Simplex convergence reached (simplex)
% -10               Optimization terminated: Stall Flights Limit reached (swarm)
% -11               Other termination status (cmaes)
% -12               Termination function change tolerance criteria reached


