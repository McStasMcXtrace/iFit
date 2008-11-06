function [pars,fval,exitflag,output] = fminhooke(fun, pars, options)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fminhooke(FUN,PARS,[OPTIONS]) Hooke-Jeeves direct search
%
% This minimization method uses Hooke-Jeeves direct search optimization
% 
% Calling:
%   fminhooke(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminhooke(fun, pars, options) same as above, with customized options (optimset)
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminhooke(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%   optimset('fminhooke')
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
% Reference: Arthur F. Kaupe Jr., Communications of the ACM, Vol 6. (1963) 313
% R. Hooke and T. A. Jeeves, Journal of the ACM, Vol. 8, April 1961, pp. 212.
% Contrib: C. T. Kelley, 1998, Iterative Methods for Optimization
%
% Version: $Revision: 1.6 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='';
  options.TolFun =1e-4;
  options.TolX   =1e-12;
  options.MaxIter='10*numberOfVariables';
  options.MaxFunEvals=1000;
  options.algorithm  = [ 'Hooke-Jeeves direct search (by Kelley) [' mfilename ']' ];
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

options=fmin_private_std_check(options, feval(mfilename,'defaults'));

if strcmp(options.Display,'iter')
  fmin_private_disp_start(mfilename, fun, pars);
end

% call the optimizer
[pars,fval,exitflag,output] = hooke(pars(:), fun, options);
output.options=options;

% private function ------------------------------------------------------------

function [pars, fval, istop, output] = hooke(x0, f, options)
%
% Hook Jeeves
%
% C. T. Kelley, July 10, 1998
%
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x, histout] = hooke(x0, f, budget, scales, tol, v)
%
% inputs:
%
%       x0= initial iterate 
%       f = objective function 
%       budget = max f evals (default=50*number of variables)
%                 The iteration will terminate after the iteration that
%                 exhausts the budget
%       scales = the decreasing sequence of stencil sizes
%                 This is an optional argument and the default is
%                 1, 1/2, ... 1/128
%       tol = termination tolerance 
%            difference is function values at successive iterations
%            (default = 1.d-6)
%       v = matrix of search directions (default $v = I$)
%       target = value of f at which the iteration should be terminated
%             This is an optional argument, which you SHOULD set to
%             something reasonable for your problem. The default is
%             as close to no limit as we can get, -1.d8
%
%
% outputs:
%
%        final results = x
%        iteration histor = histout itout x 5
%        histout = iteration history, updated after each nonlinear iteration
%                 = lhist x 2 array, the rows are
%                   [fcount, fval, sval]
%                   fcount = cumulative function evals
%                   fval = current best function value
%                   sval = current scale
%
% set debug = 1 to print iteration stats
%
debug=0;
%
%
n=length(x0); fcount=0; dir=eye(n);

tolf  = options.TolFun;
maxit = options.MaxIter;
bud   = options.MaxFunEvals;
scal  =-(0:maxit)';
scal  =2.^scal;

%
% cache control paramters
%
global cache_size cache cache_ptr cache_fvals
cache_size=4*n; cache=zeros(n,cache_size); cache_ptr=0;
cache_fvals=zeros(cache_size,1);
%
% main loop to sweep the scales
%
fcount=1; h=scal(1); fv=geval(f,x0); 
histout=[fcount,fv,h];
x=x0; nscal=length(scal); ns=0; 
pars=x;
fval=fv;
istop=0;
while (ns < nscal & fcount <= bud & ~istop)
    ns=ns+1; 
    [x, sf, lhist, fcount, fv] = hjsearch(x, f, scal(ns),dir, fcount,bud);
    histout=[histout',lhist']'; 
    pars_prev=pars;
    fval_prev=fval;
    pars=x;
    fval=fv;
    % std stopping conditions
    [istop, message] = fmin_private_std_check(pars, fval, ns, fcount, ...
      options, pars_prev);
    if strcmp(options.Display, 'iter')
      fmin_private_disp_iter(ns, fcount, f, pars, fval);
    end
%
% end main loop
%
end


% output results --------------------------------------------------------------
if istop==0, message='Algorithm terminated normally'; end
output.iterations = ns;
output.algorithm  = options.algorithm;
output.message    = message;
output.funcCount  = fcount;

if (istop & strcmp(options.Display,'notify')) | ...
   strcmp(options.Display,'final') | strcmp(options.Display,'iter')
  fmin_private_disp_final(output.algorithm, output.message, output.iterations, ...
    output.funcCount, f, pars, fval);
end

% end of main hooke function --------------------------------------------------
%
% Search with a single scale
%
function [x, sf, lhist, finc, fv] = hjsearch(xb, f, h,dir, fcount,bud)
x=xb; xc=x; sf=0; finc=fcount; lhist=[];
[x, fv, sf, numf] = hjexplore(xb, xc, f, h, dir);
finc=finc+numf;
if sf==1 thist=[finc, fv, h]; lhist=[lhist',thist']'; end
while sf==1 & finc < bud
%
% pattern move
%
    d=x-xb; xb=x; xc=x+d; fb=fv;
    [x, fv, sf, numf] = hjexplore(xb, xc, f, h, dir,fb);
    finc=finc+numf; 
    if sf == 0 % pattern move fails!
       [x, fv, sf, numf] = hjexplore(xb, xb, f, h, dir, fb);
       finc=finc+numf; 
    end
    if sf==1 thist=[finc, fv, h]; lhist=[lhist',thist']'; end
end
%
% exploratory move
%
function [x, fv, sf, numf] = hjexplore(xb, xc, f, h, dir, fbold);
global cache_size cache cache_ptr cache_fvals
n=length(xb); x=xb; numf=0;
if nargin == 5
[fb,ctr]=geval(f,x); 
numf=numf+ctr; 
else
fb=fbold; 
end
fv=fb; 
xt=xc; sf=0; dirh=h*dir;
fvold=fv;
for k=1:n
    p=xt+dirh(:,k); ft=feval(f,p); numf=numf+1;
    if(ft >= fb) p=xt-dirh(:,k); [ft,ctr]=geval(f,p); numf=numf+ctr; end
    if(ft < fb) sf=1; xt=p; fb=ft; end
end
if sf==1 x=xt; fv=fb; end
%
%
function [fs,ctr]=geval(fh,xh)
global cache_size cache cache_ptr cache_fvals
for i=1:cache_size
    nz(i)=norm2(xh-cache(:,i));
end
[vz,iz]=min(nz);
if vz == 0 & cache_ptr ~=0
    fs=cache_fvals(iz);
    ctr=0;
else
    fs=feval(fh,xh);
    ctr=1;
    cache_ptr=mod(cache_ptr,cache_size)+1;
    cache(:,cache_ptr)=xh;
    cache_fvals(cache_ptr)=fs;
end

function n=norm2(x)
x = x(:);
n=sqrt(sum(abs(x).*abs(x)));


