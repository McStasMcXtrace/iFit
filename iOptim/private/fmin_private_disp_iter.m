function fmin_private_disp_iter(iteration, funccount, fun, pars, fval)
% function called during minimization procedure
%header = ' Iteration   Func-count     min f(x)         Parameters';
if length(pars) > 5, pars=pars(1:5); end
spars=mat2str(pars(:)');  % as a row
if length(spars) > 50, spars=[ spars(1:47) ' ...' ]; end
disp(sprintf(' %5.0f        %5.0f     %12.6g         %s', iteration, funccount, fval, spars));
