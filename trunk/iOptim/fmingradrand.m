function [pars,criteria,message,output] = fmingradrand(fun, pars, options)
% [MINIMUM,FVAL,MESSAGE,OUTPUT] = FMINGRADRAND(FUN,PARS,[OPTIONS]) random gradient optimizer
%
% This optimization method uses a gradient method with random directions.
% Namely, it first determine a random direction in the optimization space
% and then uses a Newton method. This is repeated iteratively until success.
% This method is both fast and less sensitive to local minima traps.
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
%          MESSAGE return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% Contrib: Sheela V. Belur(sbelur@csc.com)
% Reference: computer Methods in Applied Mechanics & Engg, Vol  19, 99-106 , 1979
%
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='off';
  options.TolFun =1e-4;
  options.TolX   =1e-3;
  options.MaxIter=300;
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

% call the optimizer
[pars,criteria,message,output] = ossrs(fun, pars, options);

% private function ------------------------------------------------------------



function [x,fmn,message,output]=ossrs(fun, x, options)
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
if ~strcmp(options.Display,'off')
  disp('Initial Function Value and Decision Variables')
  fmn
  x
end

message='';
stdx=0.05;
nor=0;fmn0=fmn;
iterations=0; funcount=0;
while(nor<options.MaxIter)
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
  if ~strcmp(options.Display,'off')
    nor,fmn,x
  end
  
  % std stopping conditions
  [istop, message] = fmin_private_std_check(x, fmn, iterations, funcount, options);
  
  if istop
    break
  end
  
end

output.iterations = iterations;
output.algorithm  = mfilename;
output.message    = message;
output.funcCount  = funcount;


  
