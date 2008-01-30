function fmin_private_disp_final(algorithm, message, iteration, funccount, fun, pars, fval)
% function called at the end of the minimization procedure

disp([ '** Finishing minimization of ' localChar(fun) ' using algorithm ' localChar(algorithm) ]);
%display last iteration
disp(' Iteration   Func-count     min f(x)         Parameters');
fmin_private_disp_iter(-iteration, funccount, fun, pars, fval)
disp( [ '  Status: ' message ]);

