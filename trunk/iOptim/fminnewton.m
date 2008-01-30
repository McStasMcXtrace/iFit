function [pars,fval,exitflag,output] = fminnewton(fun, pars, options)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fminnewton(FUN,PARS,[OPTIONS]) Newton search
%
% This minimization method uses the Dogleg trust region, Newton model, dense algorithm 
% 
% Calling:
%   fminnewton(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminnewton(fun, pars, options) same as above, with customized options (optimset)
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminnewton(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%   optimset('fminnewton')
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
% Reference: W. Press, Numerical Recipes, Cambridge (1988)
% Contrib: C. T. Kelley, 1998, Iterative Methods for Optimization
%
% Version: $Revision: 1.1 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='off';
  options.TolFun =1e-4;
  options.TolX   =1e-6;
  options.MaxIter=20;
  options.MaxFunEvals=1000;
  pars = options;
  return
end

if nargin <= 2
  options=[];
end
if isempty(options)
  options=feval(mfilename, 'defaults');
end

options.algorithm  = [ 'Steihaug Newton-CG-Trust (by Kelley) [' mfilename ']' ];

options=fmin_private_std_check(options);

if strcmp(options.Display,'iter')
  fmin_private_disp_start(mfilename, fun, pars);
end

% call the optimizer
[pars,fval,exitflag,output] = ntrust(pars(:), fun, options);

function [pars,fval,istop,output] = ntrust(x0,f,options)
%
%
% C. T. Kelley, Dec 15, 1997
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x,histout,costdata] = ntrust(x0,f,tol,maxit,resolution)
%
% Dogleg trust region, Newton model, dense algorithm 
%
% Input: x0 = initial iterate
%        f = objective function,
%            the calling sequence for f should be
%            [fout,gout]=f(x) where fout=f(x) is a scalar
%              and gout = grad f(x) is a COLUMN vector
%        tol = termination criterion norm(grad) < tol
%        maxit = maximum iterations (optional) default = 100
%        resolution = estimated accuracy in functions/gradients (optional)
%                     default = 1.d-12
%                     The finite difference increment in the difference
%                     Hessian is set to sqrt(resolution). 
%                     
%
% Output: x = solution
%         histout = iteration history   
%             Each row of histout is
%       [norm(grad), f, TR radius, iteration count] 
%         costdata = [num f, num grad, num hess] 
%
% Requires: diffhess.m, dirdero.m
%
% set maxit, the resolution, and the difference increment
%
tol   = options.TolFun;
maxit = options.MaxIter;
resolution = 1.d-12;
hdiff=sqrt(resolution);
%
maxit=100; itc=1; xc=x0; n=length(x0);
fc = feval(f,xc);
gc = gradest(f, xc); gc=reshape(gc, size(xc));
pars=xc;
fval=fc;
numf=1; numg=1; numh=0; istop=0; message='';
% Iniitalize the TR radius, not a profound choice.

trrad=min(norm(gc),10);
%
ijob=1;
jdata=[];
while(norm(gc) > tol & itc <= maxit) & ~istop
  pars_prev=pars;
  fval_prev=fval;
  if ijob == 1
    hess=diffhess(xc,f,gc,hdiff);
    numf=numf+n; numg=numg+n; numh=numh+1;
  else
    jdata=sdata;
  end   
  itc=itc+1;
  [xp,newrad,idid,sdata,nf]=trfix(xc, f, hess, gc, fc, trrad,ijob,jdata);
  numf=numf+nf;
  ijob=idid;
  if idid == 2
     sdata=jdata;
  elseif idid == 3
     xpold=xp; trold=trrad; sdata=jdata;
  elseif idid == 4
     xp=xpold; newrad=trold; ijob=1; idid=1;
  end
  xc=xp; trrad=newrad;
  if idid==1
    fc = feval(f,xc);
    gc = gradest(f, xc); gc=reshape(gc, size(xc));
    numf=numf+1; numg=numg+1;
  end
  pars=xc; fval=fc;
	% std stopping conditions
  [istop, message] = fmin_private_std_check(pars, fval, itc, numf, ...
      options, pars_prev, fval_prev);
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

function [xp, newrad,idid,sdata,nf] =...
trfix(xc, f, hc, gc, fc, oldrad,ijob,jdata)
%
%
%     C. T. Kelley, Dec 15, 1997
%
%     This code comes with no guarantee or warranty of any kind.
%
%     function [xt, newrad, idid, sdata, nf] 
%                = trfix(xc, f, hc, gc, oldrad, ijob, jdata)
%
%     Figure out what the new trust region radius and new point are
%
%     This code is called by ntrust.m
%     There is no reason for you to call this directly.
%     
%     Input: xc = current point
%     f  = objective 
%     hc = current Hessian 
%     gc = current gradient
%     fc = current function value
%     oldrad = current TR radius
%     ijob = what to do now: 1 = fresh start
%                            2 = TR radius reduction in progress
%                            3 = attempt TR radius expansion
%     jdata = Newton direction when ijob = 1 or 2, avoids recomputation
%     nf = number of function evaluations 
%
%     Output: xp = new point
%     newrad = new TR radius
%     idid = result flag: 1 = step ok
%                         2 = TR radius reduced, step not ok
%                         3 = expansion attempted, save step and try
%                             to do better
%                         4 = expansion step failed, use the last good step
%     sdata = Newton direction to use in next call if idid > 1
%
%
%     Find the Cauchy point
%
%     bflag=1 means that the trial point is on the TR boundary and is not
%             the Newton point
%
nf=0;
bflag=0;
idid=1;
trrad=oldrad;
mu=gc'*hc*gc;
mu1=gc'*gc;
dsd=-gc; 
if ijob == 1
   dnewt=hc\dsd;
else
   dnewt=jdata;
end
sdata=dnewt;
if mu > 0
   sigma = mu1/mu;
   if(sigma*norm(gc)) > trrad
      sigma=trrad/norm(gc);
   end
   cp = xc-sigma*gc;
else
%
%     If steepest descent direction is a direction of negative curvature
%     take a flying leap to the boundary.
%
   bflag=1;
   sigma=trrad/norm(gc);
   cp=xc-sigma*gc;
end
%
%     If CP is on the TR boundary, that's the trial point.
%     If it's not, compute the Newton point and complete the dogleg.
%
if bflag==1
   xt=cp;
else
%
%     If we get to this point, CP is in the interior and the steepest
%     descent direction is a direction of positive curvature.
%
   dsd=-gc; dnewt=hc\dsd;
   xn=xc+dnewt;
   mu2=dsd'*dnewt;
%
%     If the Newton direction goes uphill, revert to CP.
%
   if mu2 <= 0
       xt=cp;
%
%     If the Newton point is inside, take it.
%
   elseif norm(dnewt) <= trrad
       xt=xn;
%
%    Newton is outside and CP is inside. Find the intersection of the
%    dog leg path with TR boundary.
%
   else
       d1=sigma*gc; d2=d1+dnewt;
       aco=d2'*d2; bco=-2*d1'*d2; cco= (d1'*d1) - trrad*trrad;
       xi=(-bco+sqrt((bco*bco) - 4*aco*cco))/(2*aco);
       xt=cp + xi*(xn-cp);
       bflag=1;
   end
end
%
%     Now adjust the TR radius using the trial point
%
st=xt-xc; ft=feval(f,xt); ared=ft-fc; nf=nf+1;
pred=gc'*st + .5* (st'*hc*st);
if ared/pred < .25
   xt=xc;
   trrad=norm(st)*.5;
   idid=2;
   if ijob == 3 idid = 4; end
elseif ared/pred > .75 & bflag==1
   trrad=trrad*2;
   idid=3;
end
newrad=trrad;
xp=xt;
