function [pars,fval,exitflag,output] = fminsimplex(fun, pars, options, constraints)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fminsimplex(FUN,PARS,[OPTIONS],[CONSTRAINTS]) Nelder-Mead Simplex
%
% This minimization method uses the Nelder-Mead Simplex
% 
% Calling:
%   fminsimplex(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminsimplex(fun, pars, options) same as above, with customized options (optimset)
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminsimplex(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%   optimset('fminsimplex')
%
%  CONSTRAINTS may be specified as a structure
%   constraints.min= vector of minimal values for parameters
%   constraints.max= vector of maximal values for parameters
%   constraints.fixed= vector having 0 where parameters are free, 1 otherwise
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
% Reference: Nelder and Mead, Computer J., 7 (1965) 308
% Contrib: C. T. Kelley, 1998, Iterative Methods for Optimization
%
% Version: $Revision: 1.9 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='';
  options.TolFun =1e-4;
  options.TolX   =1e-12;
  options.MaxIter=50;
  options.MaxFunEvals=500;
  options.algorithm  = [ 'Nelder-Mead Simplex (by Kelley) [' mfilename ']' ];
  pars = options;
  return
end

if nargin <= 2
  options=[];
end
if isempty(options)
  options=feval(mfilename, 'defaults');
end

options=fmin_private_std_check(options, feval(mfilename,'defaults'));

if strcmp(options.Display,'iter')
  fmin_private_disp_start(mfilename, fun, pars);
end

% handle constraints
if nargin <=3
  constraints=[];
elseif nargin >= 4 & isnumeric(constraints) 
  if nargin == 4,               % given as fixed index vector
    fixed = constraints; constraints=[];
    constraints.fixed = fixed;  % avoid warning for variable redefinition.
  else                          % given as lb,ub parameters (nargin==5)
    lb = constraints;
    constraints.min = lb;
    constraints.max = ub;
  end
end

if isfield(constraints, 'min')  % test if min values are valid
  index=find(isnan(constraints.min) | isinf(constraints.min));
  constraints.min(index) = -2*abs(pars(index));
  index=find(pars == 0);
  constraints.min(index) = -1;
end
if isfield(constraints, 'max')  % test if max values are valid
  index=find(isnan(constraints.max) | isinf(constraints.min));
  constraints.max(index) = 2*abs(pars(index));
  index=find(pars == 0);
  constraints.max(index) = 1;
end
if ~isfield(constraints, 'min')
  constraints.min = -2*abs(pars); % default min values
  index=find(pars == 0);
  constraints.min(index) = -1;
end
if ~isfield(constraints, 'max')
  constraints.max =  2*abs(pars); % default max values
  index=find(pars == 0);
  constraints.max(index) = 1;
end
if isfield(constraints, 'fixed') % fix some of the parameters if requested
  index = find(fixed);
  constraints.min(index) = pars(index); 
  constraints.max(index) = pars(index); % fix some parameters
end

% Initial population (random start): x) must be n x n+1
n=prod(size(pars));
prow=reshape(pars,1,n);

% Handle defaults for the step size.
Steps = constraints.max-constraints.min;
if prod(size(Steps))<n  % Only a scalar given.
    zerostep = 1e-3;
    Steps=Steps(1)*prow+zerostep*(prow==0);
end;

% The simplex is (n+1) row vectors.
x0=repmat(prow,n+1,1);
for i=1:n
    x0(i+1,i)=x0(i+1,i)+Steps(i);
end

% call the optimizer
[pars,fval,exitflag,output] = nelder(x0', fun, options);
output.options=options; output.constraints=constraints;

% private function ------------------------------------------------------------


function [pars, fval, istop, output]=nelder(x0,f,options)
%
% Nelder-Mead optimizer, No tie-breaking rule other than MATLAB's sort
%
% C. T. Kelley, December 12, 1996
%
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x,lhist,histout,simpdata] = nelder(x0,f,tol,maxit,budget)
%
% inputs:
%	vertices of initial simplex = x0 (n x n+1 matrix)
%          The code will order the vertices for you and no benefit is
%          accrued if you do it yourself.
%
%       objective function = f
%
%       termination tolerance = tol
%       maximum number of iterations = maxit (default = 100)
%           As of today, dist = | best value - worst value | < tol
%           or when maxit iterations have been taken
%       budget = max f evals (default=50*number of variables)
%                 The iteration will terminate after the iteration that
%                 exhausts the budget
%
%
% outputs:
%	final simplex = x (n x n+1) matrix
%
%       number of iterations before termination = itout (optional)
%       iteration histor = histout itout x 5
%         histout = iteration history, updated after each nonlinear iteration
%                 = lhist x 5 array, the rows are
%                   [fcount, fval, norm(grad), dist, diam]
%                   fcount = cumulative function evals
%                   fval = current best function value
%                   norm(grad) = current simplex grad norm
%                   dist = difference between worst and best values
%                   diam = max oriented length
%       simpdata = data for simplex gradient restart 
%              = [norm(grad), cond(v), bar f]
%
% initialize counters
%
lhist=0; fcount=0;

tol   = options.TolFun;
maxit = options.MaxIter;
budget= options.MaxFunEvals;
%
% set debug=1 to print out iteration stats
%
debug=0;
%
% Set the N-M parameters
%
rho=1; chi=2; gamma=.5; sigma=.5;
dsize=size(x0); n=dsize(1);
if nargin < 4 maxit=100; end
if nargin < 5 budget=100*n; end
%
% set the paramters for stagnation detection/fixup
% setting oshrink=0 gives vanilla Nelder-Mead
%
oshrink=1; restartmax=3; restarts=0;
%
%
% Order the vertices for the first time
%
x=x0; [n,m]=size(x); histout=zeros(maxit*3,5); simpdata=zeros(maxit,3);
itout=0; orth=0;
xtmp=zeros(n,n+1); z=zeros(n,n); delf=zeros(n,1);
for j=1:n+1; fv(j)=feval(f,x(:,j)); end; fcount=fcount+n+1;

[fs,is]=sort(fv); xtmp=x(:,is); x=xtmp; fv=fs;
itc=0; dist=fv(n+1)-fv(1);
diam=zeros(n,1);
for j=2:n+1
   v(:,j-1)=-x(:,1)+x(:,j);
   delf(j-1)=fv(j)-fv(1);
   diam(j-1)=norm2(v(:,j-1));
end
sgrad=v'\delf; alpha=1.d-4*max(diam)/norm2(sgrad);
lhist=lhist+1;
histout(lhist,:)=[fcount, fv(1), norm(sgrad,inf), 0, max(diam)];
%
% main N-M loop
%
istop=0;
pars=mean(x');
fval=mean(fv);
while(itc < maxit & dist > tol & restarts < restartmax & fcount <= budget & ~istop)
    fbc=sum(fv)/(n+1);
    xbc=sum(x')'/(n+1);
    sgrad=v'\delf;
    simpdata(itc+1,1)=norm2(sgrad);
    simpdata(itc+1,2)=cond(v);
    simpdata(itc+1,3)=fbc;
    if(det(v) == 0)
        istop=-11;
        message='simplex collapse';
        break
    end
    happy=0; itc=itc+1; itout=itc; how='';
%
% reflect
%
    y=x(:,1:n);
    xbart = sum(y')/n;  % centriod of better vertices
    xbar=xbart';
    xr=(1 + rho)*xbar - rho*x(:,n+1);
    fr=feval(f,xr); fcount=fcount+1;
    if(fr >= fv(1) & fr < fv(n)) happy = 1; xn=xr; fn=fr; end;
    if(happy==1) how='reflect'; end
%
% expand
%
    if(happy == 0 & fr < fv(1))
        xe = (1 + rho*chi)*xbar - rho*chi*x(:,n+1);
        fe=feval(f,xe); fcount=fcount+1;
        if(fe < fr) xn=xe;  fn=fe; happy=1; end
        if(fe >=fr) xn=xr;  fn=fr; happy=1; end
        if(happy==1) how='expand'; end
    end
%
% contract
%
   if(happy == 0 & fr >= fv(n) & fr < fv(n+1))
%
% outside contraction
%
       xc=(1 + rho*gamma)*xbar - rho*gamma*x(:,n+1);
       fc=feval(f,xc); fcount=fcount+1;
       if(fc <= fr) xn=xc; fn=fc; happy=1; end;
       if(happy==1) how='outside'; end;
   end
%
% inside contraction
%
   if(happy == 0 & fr >= fv(n+1))
       xc=(1 - gamma)*xbar+gamma*x(:,n+1);
       fc=feval(f,xc); fcount=fcount+1;
       if(fc < fv(n+1)) happy=1; xn=xc; fn=fc; end;
       if(happy==1) how='inside'; end;
   end
%
%  test for sufficient decrease, 
%  do an oriented shrink if necessary
%
   if(happy==1 & oshrink==1)
       xt=x; xt(:,n+1)=xn; ft=fv; ft(n+1)=fn;
%       xt=x; xt(:,n+1)=xn; ft=fv; ft(n+1)=feval(f,xn); fcount=fcount+1;
       fbt=sum(ft)/(n+1); delfb=fbt-fbc; armtst=alpha*norm2(sgrad)^2;
       if(delfb > -armtst/n) 
           restarts=restarts+1;
           orth=1; diams=min(diam);
           sx=.5+sign(sgrad); sx=sign(sx);
if debug==1
           [itc, delfb, armtst]
end
           happy=0;
           for j=2:n+1; x(:,j)=x(:,1); 
           x(j-1,j)=x(j-1,j)-diams*sx(j-1); end;
           how='shrink';
       end
   end
%
%  if you have accepted a new point, nuke the old point and
%  resort
%
   if(happy==1)
       x(:,n+1)=xn; fv(n+1)=fn;
%       x(:,n+1)=xn; fv(n+1)=feval(f,xn); fcount=fcount+1;
       [fs,is]=sort(fv); xtmp=x(:,is); x=xtmp; fv=fs;
   end
%
% You're in trouble now! Shrink or restart.
%
   if(restarts >= restartmax) message='stagnation in Nelder-Mead'; istop=-11; break; end
   if(happy == 0 & restarts < restartmax)
       if(orth ~=1) how='shrink'; end;
       if(orth ==1) 
       if debug == 1 disp(' restart '); end
       orth=0; end;
       for j=2:n+1;
           x(:,j)=x(:,1)+sigma*(x(:,j)-x(:,1));
           fv(j)=feval(f,x(:,j));
       end
       fcount=fcount+n;
       [fs,is]=sort(fv); xtmp=x(:,is); x=xtmp; fv=fs;
   end
%
%  compute the diameter of the new simplex and the iteration data
%
   for j=2:n+1
       v(:,j-1)=-x(:,1)+x(:,j);
       delf(j-1)=fv(j)-fv(1);
       diam(j-1)=norm2(v(:,j-1));
   end
   dist=fv(n+1)-fv(1);
   lhist=lhist+1;
   sgrad=v'\delf;
   histout(lhist,:)=[fcount, fv(1), norm(sgrad,inf), dist, max(diam)];
   pars_prev=pars;
   fval_prev=fval;
   pars=mean(x');
   fval=mean(fv);
   options.procedure=[ mfilename ': ' how ];
   % std stopping conditions
   [istop, message] = fmin_private_std_check(pars, fval, itc, fcount, ...
    options, pars_prev);
   if strcmp(options.Display, 'iter')
     fmin_private_disp_iter(itc, fcount, f, pars, fval);
   end
end

% output results --------------------------------------------------------------
if istop==0, message='Algorithm terminated normally'; end
output.iterations = itc;
output.algorithm  = options.algorithm;
output.message    = message;
output.funcCount  = fcount;

if (istop & strcmp(options.Display,'notify')) | ...
   strcmp(options.Display,'final') | strcmp(options.Display,'iter')
  fmin_private_disp_final(output.algorithm, output.message, output.iterations, ...
    output.funcCount, f, pars, fval);
end

function n=norm2(x)
x = x(:);
n=sqrt(sum(abs(x).*abs(x)));

