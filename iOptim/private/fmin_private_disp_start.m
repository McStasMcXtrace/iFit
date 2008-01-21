function fmin_private_disp_start(algorithm, fun, pars, fval)
% function called at the begining of the minimization procedure

if nargin <= 3, fval=[]; end
disp([ 'Starting minimization of ' localChar(fun) ' using algorithm ' localChar(algorithm) ]);
if length(pars) > 20, pars=pars(1:20); end
spars=mat2str(pars(:)');  % as a row
disp(  '         initial parameters:');
disp(spars);
if ~isempty(fval)
  disp(  '         initial criteria:');
  disp(fval);
end
disp(' Iteration   Func-count     min f(x)         Parameters');


