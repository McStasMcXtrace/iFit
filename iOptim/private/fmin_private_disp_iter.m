function fmin_private_disp_iter(iteration, funccount, fun, pars, fval)
% function called during minimization procedure
%header = ' Iteration   Func-count     min f(x)         Parameters';
% Displays iteration information every 5 steps, then 10 steps, then 100 steps
% or at every step if iteration is negative

if iteration > 5
  if iteration > 50 & mod(iteration,100) return;
  elseif mod(iteration,10) return; end
end
iteration=abs(iteration);
if iteration == 1
  disp(' Iteration   Func-count     min f(x)         Parameters');
end

if length(pars) > 5, pars=pars(1:5); end
spars=mat2str(pars(:)');  % as a row
if length(spars) > 50, spars=[ spars(1:47) ' ...' ]; end
disp(sprintf(' %5.0f        %5.0f     %12.6g         %s', iteration, funccount, fval, spars));
