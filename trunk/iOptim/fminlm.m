function [pars,fval,exitflag,output] = fminlm(fun, pars, options)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fminlm(FUN,PARS,[OPTIONS]) Levenberg-Maquardt search
%
% This minimization method uses the Levenberg-Maquardt steepest descent 
% in Least-Squares Sense.
% 
% Calling:
%   fminlm(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminlm(fun, pars, options) same as above, with customized options (optimset)
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminlm(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%     o=fminlm('defaults')
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
% Reference: 
% Fletcher, R., (1971) Rpt. AERE-R 6799, Harwell
% Fletcher, R., Computer Journal 1970, 13, 317-322
% Contrib: Miroslav Balda, balda AT cdm DOT cas DOT cz 2009
%
% Version: $Revision: 1.1 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display  = [];        %   no print of iterations
  options.MaxIter  = 100;       %   maximum number of iterations allowed
  options.ScaleD   = [];        %   automatic scaling by D = diag(diag(J'*J))
  options.TolFun   = 1e-7;      %   tolerace for final function value
  options.TolX     = 1e-4;      %   tolerance on difference of x-solutions
  options.MaxFunEvals=1000;
  options.algorithm  = [ 'Levenberg-Maquardt (by Balda) [' mfilename ']' ];
  options.optimizer = mfilename;
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

if strcmp(options.Display,'iter')
  fmin_private_disp_start(mfilename, fun, pars);
end

options=fmin_private_std_check(options, feval(mfilename,'defaults'));

% call the optimizer
[pars, fval, exitflag, output] = LMFsolve(fun, pars(:), options);
output.options=options;

% private function ------------------------------------------------------------

function [xf, rd, istop, output] = LMFsolve(FUN, xc, options)
% LMFSOLVE  Solve a Set of Nonlinear Equations in Least-Squares Sense.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A solution is obtained by a shortened Fletcher version of the 
% Levenberg-Maquardt algoritm for minimization of a sum of squares
%  of equation residuals. 
%
% [Xf, Ssq, CNT] = LMFsolve(FUN,Xo,Options)
% FUN     is a function handle or a function M-file name that evaluates
%         m-vector of equation residuals,
% Xo      is n-vector of initial guesses of solution,
% Options is an optional set of Name/Value pairs of control parameters 
%         of the algorithm. It may be also preset by calling:
%         Options = LMFsolve('default'), or by a set of Name/Value pairs:
%         Options = LMFsolve('Name',Value, ... ), or updating the Options
%                   set by calling
%         Options = LMFsolve(Options,'Name',Value, ...).
%
%    Name   Values {default}         Description
% 'Display'     integer     Display iteration information
%                            {0}  no display
%                             k   display initial and every k-th iteration;
% 'TolFun'      {1e-7}      norm(FUN(x),1) stopping tolerance;
% 'TolX'        {1e-7}      norm(x-xold,1) stopping tolerance;
% 'MaxIter'     {100}       Maximum number of iterations;
% 'ScaleD'                  Scale control:
%               value        D = eye(m)*value;
%               vector       D = diag(vector);
%                {[]}        D(k,k) = JJ(k,k) for JJ(k,k)>0, or
%                                   = 1 otherwise,
%                                     where JJ = J.'*J
% Not defined fields of the Options structure are filled by default values.
%
% Output Arguments:
% Xf        final solution approximation
% Ssq       sum of squares of residuals
% Cnt       >0          count of iterations
%           -MaxIter,   did not converge in MaxIter iterations

% Example:  Rosenbrock valey inside circle with unit diameter
%   R  = @(x) sqrt(x'*x)-.5;    %   A distance from the radius r=0.5
%   ros= @(x) [ 10*(x(2)-x(1)^2); 1-x(1); (R(x)>0)*R(x)*1000];
%   [x,ssq,cnt]=LMFsolve(ros,[-1.2,1],'Display',1,'MaxIter',50)
% returns   x = [0.4556; 0.2059],  ssq = 0.2966,  cnt = 18.
%
% Note:   Users with old MATLAB versions (<7), which have no anonymous
% functions implemented, should call LMFsolve with named function for
% residuals. For above example it is
%   [x,ssq,cnt]=LMFsolve('rosen',[-1.2,1]);
% where the function rosen.m is of the form
%   function r = rosen(x)
%%   Rosenbrock valey with a constraint
%   R = sqrt(x(1)^2+x(2)^2)-.5;
%%   Residuals:
%   r = [ 10*(x(2)-x(1)^2)  %   first part
%         1-x(1)            %   second part
%         (R>0)*R*1000.     %   penalty
%       ];

% Reference:
% Fletcher, R., (1971): A Modified Marquardt Subroutine for Nonlinear Least
% Squares. Rpt. AERE-R 6799, Harwell

% Miroslav Balda, 
% balda AT cdm DOT cas DOT cz
% 2007-07-02    v 1.0
% 2008-12-22    v 1.1 * Changed name of the function in LMFsolv
%                     * Removed part with wrong code for use of analytical 
%                       form for assembling of Jacobian matrix
% 2009-01-08    v 1.2 * Changed subfunction printit.m for better one, and
%                       modified its calling from inside LMFsolve.
%                     * Repaired a bug, which caused an inclination to
%                       istability, in charge of slower convergence.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   OPTIONS
    %%%%%%%
    
xf=0; rd=0; istop=0; cntfun=0;

%   LMFSOLVE(FUN,Xo,Options)
    %%%%%%%%%%%%%%%%%%%%%%%%

x  = xc(:);
lx = length(x);
FUN
r  = feval(FUN,x);             % Residuals at starting point
cntfun=cntfun+1;
%~~~~~~~~~~~~~~~~~
S  = r'*r;
epsx = options.TolX(:);
epsf = options.TolFun(:);
if length(epsx)<lx, epsx=epsx*ones(lx,1); end
J = finjac(FUN,r,x,epsx);
cntfun=cntfun+length(x);
%~~~~~~~~~~~~~~~~~~~~~~~
nfJ = 2;
A = J.'*J;                    % System matrix
v = J.'*r;

D = options.ScaleD;
if isempty(D)
    D = diag(diag(A));        % automatic scaling
    for i = 1:lx
        if D(i,i)==0, D(i,i)=1; end
    end
else
    if numel(D)>1
        D = diag(sqrt(abs(D(1:lx)))); % vector of individual scaling
    else
        D = sqrt(abs(D))*eye(lx);     % scalar of unique scaling
    end
end

Rlo = 0.25; 
Rhi = 0.75;
l=1;      lc=.75;     is=0;
cnt = 0;
ipr = options.Display;
% printit(ipr,-1);                %   Table header
d = options.TolX;               %   vector for the first cycle
maxit = options.MaxIter;        %   maximum permitted number of iterations

istop=0;

while cnt<maxit && ...          %   MAIN ITERATION CYCLE
    any(abs(d) >= epsx) && ...      %%%%%%%%%%%%%%%%%%%%
    any(abs(r) >= epsf)
    d  = (A+l*D)\v;             %   negative solution increment
    x_prev=x;
    xd = x-d;
    rd = feval(FUN,xd);
    cntfun=cntfun+1;
%   ~~~~~~~~~~~~~~~~~~~
    nfJ = nfJ+1;
    Sd = rd.'*rd;
    dS = d.'*(2*v-A*d);         %   predicted reduction

    R  = (S-Sd)/dS;
    if R>Rhi                    %   halve lambda if R too high
        l = l/2;
        if l<lc, l=0; end
    elseif R<Rlo                %   find new nu if R too low
        nu = (Sd-S)/(d.'*v)+2;
        if nu<2
            nu = 2;
        elseif nu>10
            nu = 10;
        end
        if l==0
            lc = 1/max(abs(diag(inv(A))));
            l  = lc;
            nu = nu/2;
        end
        l = nu*l;
    end
    
    cnt = cnt+1;
    %if ipr~=0 && (rem(cnt,ipr)==0 || cnt==1) %   print iteration?
    %    printit(ipr,cnt,nfJ,S,x,d,l,lc)
    %end
    % std stopping conditions
    [istop, message] = fmin_private_std_check(xd, rd, cnt, cntfun, options, x_prev);
    if strcmp(options.Display, 'iter')
      fmin_private_disp_iter(iterations, funcount, fun, x, fmn);
    end
    
    if istop
      break
    end
    
    if Sd<S
        S = Sd; 
        x = xd; 
        r = rd;
        J = finjac(FUN,r,x,epsx);
        cntfun=cntfun+length(x);
%       ~~~~~~~~~~~~~~~~~~~~~~~~~
        nfJ = nfJ+1;
        A = J'*J;       
        v = J'*r;
    end
end %   while

xf = x;                         %   final solution
rd = feval(FUN,xf);
cntfun=cntfun+1;
nfJ = nfJ+1;
Sd = rd.'*rd;
if ipr, disp(' '), end
% printit(ipr,cnt,nfJ,Sd,xf,d,l,lc)

% output results --------------------------------------------------------------
if istop==0, message='Algorithm terminated normally'; end
output.iterations = cnt;
output.algorithm  = options.algorithm;
output.message    = message;
output.funcCount  = cntfun;

if (istop & strcmp(options.Display,'notify')) | ...
   strcmp(options.Display,'final') | strcmp(options.Display,'iter')
  fmin_private_disp_final(output.algorithm, output.message, output.iterations, ...
    output.funcCount, fun, x, fmn);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FINJAC       numerical approximation to Jacobi matrix
%   %%%%%%
function J = finjac(FUN,r,x,epsx)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lx=length(x);
J=zeros(length(r),lx);
for k=1:lx
    dx=.25*epsx(k);
    xd=x;
    xd(k)=xd(k)+dx;
    rd=feval(FUN,xd);
%   ~~~~~~~~~~~~~~~~    
    J(:,k)=((rd-r)/dx);
end

