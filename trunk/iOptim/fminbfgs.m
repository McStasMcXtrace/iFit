function [pars,fval,exitflag,output] = fminbfgs(fun, pars, options)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fminbfgs(FUN,PARS,[OPTIONS]) BFGS search
%
% This minimization method uses the steepest descent/
% Broyden-Fletcher-Goldfarb-Shanno method with polynomial line search
% 
% Calling:
%   fminbfgs(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminbfgs(fun, pars, options) same as above, with customized options (optimset)
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminbfgs(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%   optimset('fminbfgs')
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
% Reference: Broyden, C. G., J. of the Inst of Math and Its Appl 1970, 6, 76-90
%   Fletcher, R., Computer Journal 1970, 13, 317-322
%   Goldfarb, D., Mathematics of Computation 1970, 24, 23-26
%   Shanno, D. F.,Mathematics of Computation 1970, 24, 647-656
% Contrib: C. T. Kelley, 1998, Iterative Methods for Optimization
%
% Version: $Revision: 1.4 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='off';
  options.TolFun =1e-4;
  options.TolX   =1e-12;
  options.MaxIter=20;
  options.MaxFunEvals=1000;
  options.algorithm  = [ 'Broyden-Fletcher-Goldfarb-Shanno (by Kelley) [' mfilename ']' ];
  pars = options;
  return
end

if nargin <= 2
  options=[];
end
if isempty(options)
  options=feval(mfilename, 'defaults');
end

if strcmp(options.Display,'iter')
  fmin_private_disp_start(mfilename, fun, pars);
end

options=fmin_private_std_check(options, feval(mfilename,'defaults'));

% call the optimizer
[pars,fval,exitflag,output] = bfgswopt(pars(:), fun, options);
output.options=options;

function [pars,fval,istop,output] = bfgswopt(x0,f,options) % tol,maxit,hess0)
%
% C. T. Kelley, July 17, 1997
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x,histout] = bfgswopt(x0,f,tol,maxit,hess0)
%
% steepest descent/bfgs with polynomial line search
% Steve Wright storage of H^-1
%
% if the BFGS update succeeds 
% backtrack on that with a polynomial line search, otherwise we use SD
%
% Input: x0 = initial iterate
%        f = objective function,
%            the calling sequence for f should be
%            [fout,gout]=f(x) where fout=f(x) is a scalar
%              and gout = grad f(x) is a COLUMN vector
%        tol = termination criterion norm(grad) < tol
%              optional, default = 1.d-6
%        maxit = maximum iterations (optional) default = 20
%         hess0 = (optional)
%            function that computes the action of the
%            initial inverse Hessian on a vector. This is optional. The
%            default is H_0 = I (ie no action). The format of hess0 is
%            h0v = hess0(v) is the action of H_0^{-1} on a vector v
%
% Output: x = solution
%         histout = iteration history   
%             Each row of histout is
%       [norm(grad), f, num step reductions, iteration count]
%         costdata = [num f, num grad, num hess] 
%                 (num hess=0 always, in here for compatibility with steep.m)
%
% At this stage all iteration parameters are hardwired in the code.

tol   = options.TolFun;
maxit = options.MaxIter;
budget= options.MaxFunEvals;

%
blow=.1; bhigh=.5;
numf=0; numg=0; numh=0;
itc=0; xc=x0;
maxarm=10; nsmax=50; debug=0;
istop=0; message='';
%
n=length(x0);
fc = feval(f,xc);
gc = gradest(f, xc); gc=reshape(gc, size(xc));
pars=xc;
fval=fc;
numf=numf+1; numg=numg+1;
go=zeros(n,1); 
alpha=zeros(nsmax,1); beta=alpha;
sstore=zeros(n,nsmax); ns=0;
%
%	dsdp = - H_c^{-1} grad_+ if ns > 0
%
while(norm(gc) > tol & itc <= maxit & ~istop)
  pars_prev=pars;
  fval_prev=fval;
	dsd=-gc;
	dsdp=-gc;
	if (ns>1)
    dsdp=bfgsw(sstore,alpha,beta,ns,dsd);
	end
%
%
% compute the direction
%
	if (ns==0) 
		dsd=-gc;
  else
	  xi=-dsdp;
	  b0=-1/(y'*s);
	  zeta=(1-1/lambda)*s+xi;
	  a1=b0*b0*(zeta'*y);
	  a1=-b0*(1 - 1/lambda)+b0*b0*y'*xi;
	  a=-(a1*s+b0*xi)'*gc;
%
%		We save go=s'*g_old just so we can use it here
%		and avoid some cancellation error
%
	  alphatmp=a1+2*a/go;
	  b=-b0*go;
%
%
	  dsd=a*s+b*xi;
  end
%
%
%
  if (dsd'*gc > -1.d-6*norm(dsd)*norm(gc))
    how='loss of descent';
	  dsd=-gc;
	  ns=0;
  end
  lambda=1; 
%
%       fixup against insanely long steps see (3.50) in the book
%
  lambda=min(1,100/(1 + norm(gc)));
  xt=xc+lambda*dsd; ft=feval(f,xt); numf=numf+1;
  itc=itc+1; 
  old=1;
  if old==0 
      goal=fc+1.d-4*(gc'*dsd); iarm=0;
      if ft > goal
               [xt,iarm,lambda]=polyline(xc,fc,gc,dsd,ft,f,maxarm);
               if iarm==-1 pars=xc; fval=fc;
                 message='line search failure'; istop=-11; break
               end
      end
  end
  if old==1
   	iarm=0; goalval=.0001*(dsd'*gc);
   	q0=fc; qp0=gc'*dsd; lamc=lambda; qc=ft;
    while(ft > fc + lambda*goalval )
	    iarm=iarm+1;
      if iarm==1
         lambda=polymod(q0, qp0, lamc, qc, blow, bhigh);
      else
         lambda=polymod(q0, qp0, lamc, qc, blow, bhigh, lamm, qm);
      end
      qm=qc; lamm=lamc; lamc=lambda;
	    xt=xc+lambda*dsd;
	    ft=feval(f,xt); qc=ft; numf=numf+1;
	    if(iarm > maxarm) 
        pars=xc; fval=fc;
	      message='too many backtracks in BFGS line search'; 
	      istop=-11; break; end
    end
  end
  if istop, break; end
	s=xt-xc; y=gc; go=s'*gc;
%        lambda=norm(s)/norm(dsd);
	xc=xt; 
  fc = feval(f,xc);
  gc = gradest(f, xc); gc=reshape(gc, size(xc));
  y = gc-y; numf=numf+1; numg=numg+1;
%
%   restart if y'*s is not positive or we're out of room
%
	if (y'*s <= 0) | (ns==nsmax) 
    message='loss of positivity or storage'; 
		ns=0;
	else
		ns=ns+1; sstore(:,ns)=s;
		if(ns>1)
			alpha(ns-1)=alphatmp;
			beta(ns-1)=b0/(b*lambda);
		end
	end
	pars=xc; fval=fc;
	% std stopping conditions
	options.procedure=message;
  [istop, message] = fmin_private_std_check(pars, fval, itc, numf, ...
      options, pars_prev);
  if strcmp(options.Display, 'iter')
    fmin_private_disp_iter(itc, numf, f, pars, fval);
  end
end

% output results --------------------------------------------------------------
if istop==0, message='Algorithm terminated normally'; end
output.iterations = itc;
output.algorithm  = options.algorithm;
output.message    = message;
output.funcCount  = numf;

if (istop & strcmp(options.Display,'notify')) | ...
   strcmp(options.Display,'final') | strcmp(options.Display,'iter')
  fmin_private_disp_final(output.algorithm, output.message, output.iterations, ...
    output.funcCount, f, pars, fval);
end

%
% bfgsw
%
% C. T. Kelley, Dec 20, 1996
%
% This code comes with no guarantee or warranty of any kind.
%
% This code is used in bfgswopt.m 
% 
% There is no reason to ever call this directly.
%
% form the product of the bfgs approximate inverse Hessian
% with a vector using the Steve Wright method
%
function dnewt=bfgsw(sstore,alpha,beta,ns,dsd,hess0)
dnewt=dsd; 
if (ns<=1) return; end;
dnewt=dsd; n=length(dsd);
sigma=sstore(:,1:ns-1)'*dsd; gamma1=alpha(1:ns-1).*sigma;
gamma2=beta(1:ns-1).*sigma;
gamma3=gamma1+beta(1:ns-1).*(sstore(:,2:ns)'*dsd);
delta=gamma2(1:ns-2)+gamma3(2:ns-1);
dnewt=dnewt+gamma3(1)*sstore(:,1)+gamma2(ns-1)*sstore(:,ns);
if(ns <=2) return; end
dnewt=dnewt+sstore(1:n,2:ns-1)*delta(1:ns-2);
%
