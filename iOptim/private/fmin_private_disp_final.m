function fmin_private_disp_final(algorithm, message, iteration, funccount, fun, pars, fval)
% function called at the end of the minimization procedure

disp([ 'Finishing minimization of ' localChar(fun) ' using algorithm ' localChar(algorithm) ]);
disp([ 'Required ' num2str(iteration) ' iterations and ' num2str(funccount) ' function evaluations.'])
if length(pars) > 20, pars=pars(1:20); end
spars=mat2str(pars(:)');  % as a row
disp(  '         final parameters:');
disp(spars);
disp(  '         final criteria:');
disp(fval);
disp(  '         final status:');
disp(message);



