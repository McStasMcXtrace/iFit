function [istop, message] = fmin_private_std_check(pars, fval, iterations, funccount, options)


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
  
  if options.FunValCheck & any(isnan(fval) | isinf(fval))
    istop=-4;
    message = 'Function value is Inf or Nan (options.FunValCheck)';
  end

  if ~isempty(options.OutputFcn)
    optimValues = options;
    optimValues.iterations = iterations;
    optimValues.funccount  = funccount;
    optimValues.fval       = fval;
    istop = feval(options.OutputFcn, pars, optimValues, 'iter');
    if istop, 
      istop=-5;
      message = 'Algorithm was terminated by the output function (options.OutputFcn)';
    end
  end
