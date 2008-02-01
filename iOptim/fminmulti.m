function [pars,fval,exitflag,output] = fminmulti(fun, pars, options, constraints)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fminmulti(FUN,PARS,[OPTIONS]) Multi directional search
%
% This minimization method uses the Multi directional search, which is not very
% accurate but may be used prior to other time coonsuming techniques, to guess
% starting parameters.
% 
% Calling:
%   fminmulti(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminmulti(fun, pars, options) same as above, with customized options (optimset)
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminmulti(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%   optimset('fminmulti')
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
% Version: $Revision: 1.3 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='off';
  options.TolFun =1e-6;
  options.TolX   =1e-12;
  options.MaxIter=400;
  options.MaxFunEvals=2000;
  options.algorithm  = [ 'Multidirectional search (by Kelley) [' mfilename ']' ];
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
  constraints.min(index) = 1;
end
if ~isfield(constraints, 'min')
  constraints.min = -2*abs(pars); % default min values
  index=find(pars == 0);
  constraints.min(index) = -1;
end
if ~isfield(constraints, 'max')
  constraints.max =  2*abs(pars); % default max values
  index=find(pars == 0);
  constraints.min(index) = 1;
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
[pars,fval,exitflag,output] = mds(x0', fun, options);

% private function ------------------------------------------------------------
function [pars,fval,istop,output]=mds(x0,f,options)
%
% Multidirectional search
%
% C. T. Kelley, July 17, 1998
%
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x,lhist,histout] = mds(x0,f,tol,maxit,budget)
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
%       lhist = number of iterations before termination
%       iteration histor = histout itout x 4
%         histout = iteration history, updated after each nonlinear iteration
%                 = lhist x 4 array, the rows are
%                   [fcount, fval, norm(grad), dist, diam]
%                   fcount = cumulative function evals
%                   fval = current best function value
%                   dist = difference between worst and best values
%                   diam = max oriented length
%
% initialize counters
%
lhist=0; fcount=0;
%
% set debug=1 to print out iteration stats
%
debug=0;
tol   = options.TolFun;
maxit = options.MaxIter;
budget= options.MaxFunEvals;
%
% Set the MDS parameters
%
de=2; dc=.5; 
%
% cache control paramters
%
global cache_size cache cache_ptr cache_fvals
[n,m]=size(x0); 
cache_size=4*n; cache=zeros(n,cache_size); cache_ptr=0;
cache_fvals=zeros(cache_size,1);
%
% Order the vertices for the first time
%
x=x0; 
histout=[];
itout=0; orth=0; fv=zeros(n+1,1);
xtmp=zeros(n,n+1); z=zeros(n,n); delf=zeros(n,1);
for j=1:n+1; fv(j)=geval(f,x(:,j)); 
end; fcount=fcount+n+1;
[fs,is]=sort(fv); xtmp=x(:,is); x=xtmp; fv=fs;
itc=0; dist=fv(n+1)-fv(1);
diam=zeros(n,1);
for j=2:n+1
   v(:,j-1)=-x(:,1)+x(:,j);
   diam(j-1)=norm(v(:,j-1));
end
lhist=lhist+1;
thist=[fcount, fv(1), dist, max(diam)]; histout=[histout',thist']';
%
% main MDS loop
%
istop=0;
pars=mean(x');
fval=mean(fv);
while(itc < maxit & fval > tol & fcount <= budget & ~istop)
  happy=0;
  itc=itc+1;
  %
  %  reflect
  %
  for j=2:n+1 xr(:,j)=x(:,1) - v(:,j-1); 
    [fr(j),ctr] = geval(f,xr(:,j)); 
    fcount=fcount+ctr;
  end
  how='reflect';
  %fcount=fcount+n; 
  fvr=min(fr(2:n+1)); 
  if fv(1) > fvr
    happy=1; how='expand';
%
% expand
%
    for j=2:n+1 xe(:,j)=x(:,1) - de*v(:,j-1); fe(j) = feval(f,xe(:,j)); end
    fcount=fcount+n; fve=min(fe(2:n+1));
    if fvr > fve
        for j=2:n+1 x(:,j)=xe(:,j); fv(j)=fe(j); end
    else
        for j=2:n+1 x(:,j)=xr(:,j); fv(j)=fr(j); end
    end
  end
%
% contract
%
  if happy==0
    for j=2:n+1 x(:,j)=x(:,1) + dc*v(:,j-1); 
    [fv(j),ctr] = geval(f,x(:,j)); fcount=fcount+ctr;
    end
    how='contract';
%   fcount=fcount+n;
  end
%
%  sort the new vertices
%
  [fs,is]=sort(fv); xtmp=x(:,is); x=xtmp; fv=fs;
%
%  compute the diameter of the new simplex and the iteration data
%
  for j=2:n+1
    v(:,j-1)=-x(:,1)+x(:,j);
    diam(j-1)=norm(v(:,j-1));
  end
  dist=fv(n+1)-fv(1);
  lhist=lhist+1;
  thist=[fcount, fv(1), dist, max(diam)]; histout=[histout',thist']';

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
end % while

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
%
%
%
function [fs,ctr]=geval(fh,xh)
global cache_size cache cache_ptr cache_fvals
for i=1:cache_size
    nz(i)=norm(xh-cache(:,i));
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
