function [pars,fval,exitflag,output] = fmingradrand(fun, pars, options)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = FMINGRADRAND(FUN,PARS,[OPTIONS]) random gradient optimizer
%
% This minimization method uses a gradient method with random directions.
% Namely, it first determine a random direction in the optimization space
% and then uses a Newton method. This is repeated iteratively until success.
% This method is both fast and less sensitive to local minima traps.
% 
% Calling:
%   fmingradrand(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fmingradrand(fun, pars, options) same as above, with customized options (optimset)
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fmingradrand(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%   optimset('fmingradrand')
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% Reference: Computer Methods in Applied Mechanics & Engg, Vol  19, (1979) 99
% Contrib: Sheela V. Belur(sbelur@csc.com) 1998
%
% Version: $Revision: 1.10 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='';
  options.TolFun =1e-4;
  options.TolX   =1e-12;
  options.MaxIter=300;
  options.MaxFunEvals=1000;
  options.algorithm  = [ 'Random Gradient (by Belur) [' mfilename ']' ];
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

% call the optimizer
[pars,fval,exitflag,output] = ossrs(fun, pars, options);
output.options=options;

% private function ------------------------------------------------------------

function [x,fmn,istop,output]=ossrs(fun, x, options)
%OSSRS: Optimized Step Size Random Search
% Author of method and code: Sheela V. Belur(sbelur@csc.com)
% Reference: computer Methods in Applied Mechanics & Engg, Vol  19, 99-106 , 1979
%
% function [x,fmn]=ossrs(x)
%
% n is the number of decision variables.
% x is the starting values and returns the optimal values
% fmn is the optimal function value 
n=length(x);
f(1)=feval(fun,x);
fmn=f(1);
if strcmp(options.Display,'iter')
  fmin_private_disp_start(mfilename, fun, x, fmn);
end

istop=0;
stdx=sqrt(options.TolFun);
nor=0;fmn0=fmn;
iterations=0; funcount=0;
while(nor<options.MaxIter)
  x_prev=x;
  f_prev=fmn;
  iterations=iterations+1;
  r=randn(size(x))*stdx;
  an2=x-r;
  an3=x+r;
  f(2)=feval(fun,an2);
  f(3)=feval(fun,an3);
  nor=nor+2;
  funcount=funcount+2;
  a=f(2)-2.0*f(1)+f(3);
  if(a>0.0)
    amda=-0.5*(f(3)-f(2))/a;
    an4=x+amda*r;
    f(4)=feval(fun,an4);
    nor=nor+1;
    funcount=funcount+1;
  end

  [fmn,mn]=min(f);
  if(mn==2)    x=an2;
  elseif(mn==3)x=an3;
  elseif(mn==4)x=an4;
  end
  f(1)=fmn;
  difff=abs(fmn-fmn0);
  if(difff==0.0) nor=nor+1;
  else 
    nor=0;fmn0=fmn;
  end
  
  % std stopping conditions
  [istop, message] = fmin_private_std_check(x, fmn, iterations, funcount, options, x_prev);
  if strcmp(options.Display, 'iter')
    fmin_private_disp_iter(iterations, funcount, fun, x, fmn);
  end
  
  if istop
    break
  end
  
end

% output results --------------------------------------------------------------
if istop==0, message='Algorithm terminated normally'; end
output.iterations = iterations;
output.algorithm  = options.algorithm;
output.message    = message;
output.funcCount  = funcount;

if (istop & strcmp(options.Display,'notify')) | ...
   strcmp(options.Display,'final') | strcmp(options.Display,'iter')
  fmin_private_disp_final(output.algorithm, output.message, output.iterations, ...
    output.funcCount, fun, x, fmn);
end

