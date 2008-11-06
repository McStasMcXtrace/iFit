function [pars,fval,exitflag,output] = fminimfil(fun, pars, options)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fminimfil(FUN,PARS,[OPTIONS]) Unconstrained Implicit filtering 
%
% Implicit filtering solves unconstrained optimization problems
% Minimization of noisy functions
%
% Calling:
%   fminimfil(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminimfil(fun, pars, options) same as above, with customized options (optimset)
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminimfil(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%   fminimfil('defaults')
%  options.Hybrid specifies the algorithm to use for local hybrid optimizations.
%   This is a string with possible values 'sr1','bfgs','none'.
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% Reference: C. T. Kelley, Iterative Methods for Optimization, no. 18 in 
%   Frontiers in Applied Mathematics, SIAM, Philadelphia, 1999.
% Contrib: C. T. Kelley, 1998, Iterative Methods for Optimization
%
% Version: $Revision: 1.6 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset;
  % add Matlab std options.
  options.Display='';
  options.TolFun =1e-6;
  options.TolX   =1e-12;
  options.MaxIter=100;
  options.MaxFunEvals=1000;
  options.Hybrid = 'BFGS';
  options.algorithm  = [ 'Unconstrained Implicit filtering (by Kelley) [' mfilename ']' ];
  pars = options;
  return
end

if nargin <= 2
	options=[];
end
if isempty(options)
  options=feval(mfilename, 'defaults');
end
n = prod(size(pars));
numberOfVariables = n;
if ischar(options.MaxFunEvals), 
  options.MaxFunEvals = eval(options.MaxFunEvals); 
end

if ischar(options.MaxIter), 
  options.MaxIter = eval(options.MaxIter); 
end
if ~isfield(options,'Hybrid'), options.Hybrid=''; end
if isempty(options.Hybrid),    options.Hybrid='BFGS'; end

if strcmp(options.Display,'iter')
  fmin_private_disp_start(mfilename, fun, pars);
end

options=fmin_private_std_check(options, feval(mfilename,'defaults'));
options.algorithm  = [ 'Unconstrained Implicit filtering (by Kelley) [' mfilename '/' options.Hybrid ']' ];

% call the optimizer
[pars,fval,exitflag,output] = imfil(pars(:), fun, options);
output.options=options;

% PRIVATE original code -------------------------------------------------------
function [pars,fval,istop,output] = imfil(x0,f,options)
%
%
% C. T. Kelley, January 9, 1998
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x,fcount,histout] = imfil(x0,f,budget,scales,parms)
%
% Unconstrained implicit filtering code
% 
% IMPLICIT FILTERING with SR1 and BFGS quasi-Newton methods
%
% Input: x0 = initial iterate
%        f = objective function,
%            the calling sequence for f should be [fout]=f(x)
%        budget = max f evals 
%                 The iteration will terminate after the iteration that
%                 exhausts the budget, default=50*number of variables
%        scales = the decreasing sequence of difference increments 
%                 This is an optional argument and the default is
%                 1, 1/2, ... 1/128
%        parms = optional argument = array of conrol praamters
%
%        parms(1) = 
%             target = value of f at which the iteration should be terminated
%                 This is an optional argument, which you SHOULD set to
%                 something reasonable for your problem. The default is
%                 as close to no limit as we can get, -1.d8
%
%        parms(2) = 0 for centered diffs, 1 for forward diffs
%                   default and recommended value = 0
% 
%        parms(3) = quasi-Newton method selection
%                   0 = none, 1 = bfgs, 2 = SR1
%                   default and recommend value = 1
%
%
% Output: x = estimated minimizer
%         lhist = number of nonzero rows in histout 
%               = number of nonlinear iterations for all the scales        
%         histout = iteration history, updated after each nonlinear iteration 
%                 = lhist x 5 array, the rows are
%                   [fcount, fval, norm(sgrad), norm(step), iarm]
%                   fcount = cumulative function evals
%                   fval = current function value
%                   norm(sgrad) = current simplex grad norm
%                   norm(step) = norm of last step 
%                   iarm=line searches to date
%                        =-1 means first iterate at a new scale 
%
% This code uses centered difference approximations to the gradient. Setting
% fdiff = 1 will change this to forward differences. We do not recommend 
% that.
%
% set debug = 1 to print iteration stats
%
debug=0;
%
fcount=0; 
%
% And now for the knobs. Implicit filtering has too many of these and they
% can make a difference. The ones we use are:
%
% min_gscal (default = .01)
%   if norm(difference_grad) < min_gscal*h  we terminate at the scale
%   with success
%
% maxit and maxitarm (defaults 2 and 5)
%    At most maxit*n iterations are taken for each scale and at most
%    maxitarm step length reductions are allowed
%
% nterm controls termination on stencil failure for centered diffs (defalt = 0)
%       = 0 to terminate on stencil failure before starting the line search
%       = 1 to ignore stencil failure
%
% iquit (default = 3)
%    After iquit consecutive line search failures, we terminate the iteration
%
% beta (default = .5) 
%    step size reduction factor in the line search
%
% set the knobs
%
% min_gscal=.5; 
min_gscal=.01; 
% maxit=10; maxitarm=10; iquit=3; beta=.1;
maxit=options.MaxIter; 
maxitarm=10; iquit=3; beta=.5; nterm=0;

switch lower(options.Hybrid)
case 'sr1'
  quasi=2;
case 'none'
  quasi=0;
case 'bfgs'
  quasi=1;
otherwise % BFGS
  quasi=1;
end

flim=options.MaxFunEvals;
%
% set up the difference scales
%
quasi=1; fdiff=0;
ftol=options.TolFun;
nterm=fdiff+nterm;

dscal=-(0:maxit)'; dscal=2.^dscal;

nscal=length(dscal);
n=length(x0);
%
% sweep through the scales
%
x=x0; xold=x0; n=length(x0); v=eye(n); xc=x0; hess=eye(n); ns=0; iquitc=0;
iterations=0;
% for ns=1:nscal
while (ns < nscal & fcount <= flim & iquitc < iquit)
  ns=ns+1;
  itc=0; h=dscal(ns); z0=x; fval=feval(f,x); fcount=fcount+1;
  pars_prev=x;
  fval_prev=fval;
  stol=min_gscal*h; iarm=0; lok=1;
  [sgrad,fb,xb,sflag] =simpgrad(x,f,h*v,fval,fdiff);
  fcount=fcount+(2-fdiff)*n; 
  if norm(sgrad,inf) < stol | (sflag+nterm)==0
%
%   Convergence at this scale on stencil failure or tolerance match
%
    gc=sgrad;
    if (sflag+nterm) ~= 0 
       iquitc=iquitc+1; 
    else
       iquitc=0;
    end
    iterations=iterations+1;
  else
%
%   Take a few quasi-Newton iterates
%
    iquitc=0;
    while itc < maxit*n & fval > ftol & norm(sgrad,inf) >= stol...
        &lok==1 & fcount < flim & sflag+nterm > 0
      itc=itc+1;
      iterations=iterations+1;
%
%     compute the difference gradient, scale it
%
      gc=sgrad; 
      if(itc > 1)
        [sgrad,fb,xb,sflag]=simpgrad(x,f,h*v,fval,fdiff); 
        fcount=fcount+(2-fdiff)*n; 
      end
      dgrad=sgrad;
%
%     watch out for stencil failure!
%
      if sflag+nterm > 0
%
%     update iterate and Hessian 
%
        if itc > 1 & quasi > 0
          if quasi==1
             hess = bfupdate(x, xc, sgrad, gc, hess); 
          else
            hess = sr1up(x, xc, sgrad, gc, hess); 
          end
        end; 
        xc=x;
%
%     new direction and line search
%
        if quasi > 0
          sdir=hess\dgrad;
        else
          sdir=dgrad;
        end
        [fct, x, fval, hess, iarm]=...
              linearm(f, sdir, fval, x, hess, maxitarm,beta,h,quasi,fdiff);
        fcount=fcount+fct;
%
%     reduce scale upon failure of line search
%
        if iarm >= maxitarm
           lok=0;
           x=xb; fval=fb;
        end
      end
%
%      keep the records
%
      stepn=norm(xold-x,inf); xold=x;
%
    end % end of nonlinear step
%
  end % end of sweep through the scale
  pars=x;
  % std stopping conditions
  [istop, message] = fmin_private_std_check(pars, fval, iterations, fcount, ...
    options, pars_prev);
  if strcmp(options.Display, 'iter')
    fmin_private_disp_iter(iterations, fcount, f, pars, fval);
  end
  if istop
    break
  end
end % end of while loop over the scales

if iquitc >= iquit
  istop=-11; message='line search failures';
end

% output results --------------------------------------------------------------
if istop==0, message='Algorithm terminated normally'; end
output.iterations = iterations;
output.algorithm  = options.algorithm;
output.message    = message;
output.funcCount  = fcount;

if (istop & strcmp(options.Display,'notify')) | ...
   strcmp(options.Display,'final') | strcmp(options.Display,'iter')
  fmin_private_disp_final(output.algorithm, output.message, output.iterations, ...
    output.funcCount, f, pars, fval);
end
% end main imfil

%
%   BFGS update of Hessian; nothing fancy
%
function hess = bfupdate(x, xc, sgrad, gc, hess)
y=sgrad-gc; s=x-xc; z=hess*s;
if y'*s > 0
   hess = hess + (y*y'/(y'*s)) - (z*z'/(s'*z));
end
%
% SR1 update
%
function hess = sr1up(x, xc, sgrad, gc, hess)
y=sgrad-gc; s=x-xc; z=y - hess*s;
if z'*s ~=0
	ptst=z'*(hess*z)+(z'*z)*(z'*z)/(z'*s); 
	if ptst > 0 hess = hess + (z*z')/(z'*s); end
end
%
%    Line search for implicit filtering
%
function [fct, x, fval, hessp, iarm]=...
                linearm(f, sdir, fold, xc, hess, maxitarm,beta,h,quasi,fdiff)
lambda=1;
n=length(xc);
hessp=hess;
iarm=-1;
fct=0;
aflag=1;
dd=sdir;
smax=10*min(h,1); if norm(dd) > smax dd=smax*dd/norm(dd); end
x=xc;
fval=fold;
while iarm < maxitarm & aflag==1
    d=-lambda*dd; 
    iarm=iarm+1;
    xt=x+d; ft=feval(f,xt); fct=fct+1;
    if ft < fval & aflag==1; aflag=0; fval=ft; x=xt; end
    if aflag==1; lambda=beta*lambda; end
end
if iarm == maxitarm & aflag == 1
       % disp(' line search failure'); [iarm, h, quasi, fdiff]
%      hessp=eye(n);
end
