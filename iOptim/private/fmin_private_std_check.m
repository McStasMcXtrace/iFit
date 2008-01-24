function [istop, message] = fmin_private_std_check(pars, fval, iterations, funccount, options, pars_prev)

  istop=0; message='';

  if options.TolFun & fval <= options.TolFun
    istop=-1;
    message = [ 'Termination function tolerance criteria reached (options.TolFun=' ...
              num2str(options.TolFun) ')' ];
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
  
  if nargin >= 6
    if all(abs(pars(:)-pars_prev(:)) < abs(options.TolX*pars(:)))
      istop=-5;
      message = [ 'Termination parameter tolerance criteria reached (options.TolX=' ...
            num2str(options.TolX) ')' ];
    end
  end

  if ~isempty(options.OutputFcn)
    optimValues = options;
    if ~isfield(optimValues,'state')
      if istop,               optimValues.state='done';
      elseif iterations <= 1, optimValues.state='init';
      else                    optimValues.state='iter'; end
    end
    optimValues.iterations = iterations;
    optimValues.funccount  = funccount;
    optimValues.fval       = fval;
    istop = feval(options.OutputFcn, pars, optimValues, optimValues.state);
    if istop, 
      istop=-6;
      message = 'Algorithm was terminated by the output function (options.OutputFcn)';
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


