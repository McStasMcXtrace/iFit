function [istop, message] = fmin_private_std_check(pars, fval, iterations, funccount, options, pars_prev)
% standard checks
% fmin_private_std_check(pars, fval, iterations, funccount, options
% or
% fmin_private_std_check(pars, fval, iterations, funccount, options, pars_prev)
% or
% options=fmin_private_std_check(options);
% or
% options=fmin_private_std_check(options, default_options);


  istop=0; message='';
  
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

  % normal terminations: tolerance reached
  if options.TolFun & fval <= options.TolFun
    if nargin >= 6
      if options.TolX & all(abs(pars(:)-pars_prev(:)) < abs(options.TolX*pars(:)))
        istop=-5;
        message = [ 'Termination parameter tolerance criteria reached (options.TolX=' ...
              num2str(options.TolX) ', options.TolFun=' ...
              num2str(options.TolFun) ')' ];
      end
    end
    if ~istop
      istop=-1;
      message = [ 'Termination function tolerance criteria reached (options.TolFun=' ...
                num2str(options.TolFun) ')' ];
    end
  end
  
  % abnormal terminations
  % last check before returning a failed state
  if ~istop
    if (iterations >= floor(options.MaxIter*0.95)  | funccount >= floor(options.MaxFunEvals*.95) ) & ...
    options.TolFun & fval <= options.TolFun*2
      istop=-1;
      message = [ 'Termination function tolerance criteria nearly reached (options.TolFun=' ...
                num2str(options.TolFun*2) ];
    end

    if options.MaxIter & iterations >= options.MaxIter
      istop=-2;
      message = [ 'Maximum number of iterations reached (options.MaxIter=' ...
                num2str(options.MaxIter) ')' ];
    end

    if options.MaxFunEvals & funccount >= options.MaxFunEvals
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


