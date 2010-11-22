function fmin_private_disp_iter(options, iteration, funccount, fun, pars, fval)
% function called during minimization procedure
%
% Displays iteration information every 5 steps, then 10 steps, then 100 steps
% or at every step if iteration is negative


if iteration > 5
  if iteration > 50 & mod(iteration,100) return;
  elseif mod(iteration,10) return; end
end
iteration=abs(iteration);

if isfield(options,'Display')
  if (strcmp(options.Display, 'iter') || strcmp(options.Display, 'start')) & iteration == 1
    disp(' Iteration   Func_count     min[f(x)]        Parameters');
    options.Display='iter'; % display first line
  elseif strcmp(options.Display, 'iter')
    if length(pars) > 20, pars=pars(1:20); end
    spars=mat2str(pars(:)', 4);  % as a row
    if length(spars) > 50, spars=[ spars(1:47) ' ...' ]; end
    disp(sprintf(' %5.0f        %5.0f     %12.6g         %s', iteration, funccount, fval, spars));
  end
end

