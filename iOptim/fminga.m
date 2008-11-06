function [pars,fval,exitflag,output] = fminga(fun, pars, options, constraints, ub)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = FMINGA(FUN,PARS,[OPTIONS],[CONSTRAINTS]) genetic algorithm optimizer
%
% This minimization method uses a Genetic Algorithm (real coding), optionally 
% with constraints on function parameters.
% 
% Calling:
%   fminga(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminga(fun, pars, options) same as above, with customized options (optimset)
%   fminga(fun, pars, options, fixed) 
%     is used to fix some of the parameters. The 'fixed' vector is then 0 for
%     free parameters, and 1 otherwise.
%   fminga(fun, pars, options, lb, ub) 
%     is used to set the minimal and maximal parameter bounds, as vectors.
%   fminga(fun, pars, options, constraints) 
%     where constraints is a structure (see below).
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminga(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%   optimset('fminga')
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
%
% Contrib:
% By: Javad Ivakpour javad7@gmail.com, May 2006
%
% Version: $Revision: 1.12 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='';
  options.TolFun =1e-4;
  options.TolX   =1e-12;
  options.MaxIter=1000;
  options.MaxFunEvals=1000*100;
  options.PopulationSize=50;
  options.algorithm  = [ 'Genetic Algorithm (real coding by Ivakpour) [' mfilename ']' ];
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

options=fmin_private_std_check(options, feval(mfilename,'defaults'));

% call the optimizer
[pars,fval,exitflag,output] = GA(fun, pars(:)', options, constraints);
output.options=options; output.constraints=constraints;

% private function ------------------------------------------------------------


function [pars, fval, istop, output] = GA(fun, pars, options, constraints)
% Genetic Algorithm(real coding)
% By: Javad Ivakpour
% E-mail: javad7@gmail.com
% May 2006
% Goal: find maximum of function that introduced in fun00.m file in current
% directory and can be plotted in plot00
% this file is also include the random search for comparision

%------------------------        parameters        ------------------------
% befor using this function you must specified your function in fun00.m
% file in current directory and then set the parameters
var=length(pars);            % Number of variables (this item must be equal to the
                  %   number of variables that is used in the function
n=options.PopulationSize;            % Number of population

nmutationG=20;                  %number of mutation children(Gaussian)
nmutationR=20;                  %number of mutation children(random)
nelit=2;                        %number of elitism children
valuemin=constraints.min(:);       % min possible value of variables
valuemax=constraints.max(:);       % max possible value of variables
pars=pars(:);
fval=feval(fun, pars);
%-------------------------------------------------------------------------
nmutation= nmutationG+nmutationR;
max1     = zeros(nelit,var);
parent   = zeros(n,var);
p        = zeros(n,var);
sigma=(valuemax-valuemin)/10;    %Parameter that related to Gaussian
                                 %   function and used in mutation step
max1=zeros(nelit,var);
parent=zeros(n,var);
for l=1:var
    p(:,l)=valuemin(l)+rand(n,1).*(valuemax(l)-valuemin(l));
end
p(1,:)=pars;  % the starting configuration is the first guess
g=0;
istop=0;
iterations=0;
funcount  =0;
istop     =0;

if strcmp(options.Display,'iter')
  fmin_private_disp_start(mfilename, fun, pars);
end

%-------------   ****    termination criteria   ****-------------
while 1
    sigma=sigma./(1.05);% reducing the sigma value
    %  ------  **** % reducing the number of mutation()random   **** ----
    g=g+1;
    if g>10 & nmutationR>0
        g=0;
        nmutationR=nmutationR-1;
        nmutation=nmutationG+nmutationR;
    end


    %-------------   ****    function evaluation   ****-------------
    for i=1:n
        y(i)=-feval(fun,p(i,:));  % search maximum of -function
    end
    maxy=max(y);
    miny=min(y);
    
    s=sort(y);
    fvalsorted(1:nelit)=s(n:-1:n-nelit+1);
    if nelit==0
        fvalsorted(1)=s(n);
        for i=1:n
            if y(i)==fvalsorted(1)
                max1(1,:)=p(i,:);
            end
        end
    end
    for k=1:nelit
        for i=1:n
            if y(i)==fvalsorted(k)
                max1(k,:)=p(i,:);
            end
        end
    end
    
    y=y-miny*1.02;
    sumd=y./sum(y);


    %-------------   ****   Selection: Roulette wheel   ****-------------
    for l=1:n
        sel=rand;
        sumds=0;
        j=1;
        while sumds<sel
            sumds=sumds+sumd(j);
            j=j+1;
        end
        parent(l,:)=p(j-1,:);
    end
    p=zeros(n,var);

    %-------------   ****    regeneration   ****-------------
    for l=1:var


        %-------------   ****    cross-over   ****-------------
        for j=1:ceil((n-nmutation-nelit)/2)
            t=rand*1.5-0.25;
            p(2*j-1,l)=t*parent(2*j-1,l)+(1-t)*parent(2*j,l);
            p(2*j,l)=t*parent(2*j,l)+(1-t)*parent(2*j-1,l);
        end


        %-------------   ****    elitism   ****-------------
        for k=1:nelit
            p((n-nmutation-k+1),l)=max1(k,l);
        end


        %-------------   ****    mutation (Gaussian)  ****-------------
        phi=1-2*rand(nmutation-nmutationR,1);
        z = erfinv(phi)*sqrt(2);
        p(n-nmutation+1:n-nmutationR,l) = z * sigma(l)+parent(n-nmutation+1:n-nmutationR,l);
        
        %-------------   ****    mutation (Random)  ****-------------
        for i=n-nmutationR+1:n
            p(i,1:var)=(valuemin(1:var)+rand(var,1).*(valuemax(1:var)...
                -valuemin(1:var)))';
        end
        for i=1:n % constraints parameters within range
          if p(i,l)<valuemin(l),    p(i,l)=valuemin(l);
          elseif p(i,l)>valuemax(l),p(i,l)=valuemax(l);
          end
        end
    end
    
    % std stopping conditions
    fval_prev = fval;
    fval = -fvalsorted(1);
    pars_prev=pars(:)';
    pars =  max1(1,:); pars = pars(:)';
    iterations = iterations+1;
    funcount = n*iterations;
    [istop, message] = fmin_private_std_check(pars, fval, iterations, funcount, options, pars_prev);
    if strcmp(options.Display, 'iter')
      fmin_private_disp_iter(iterations, funcount, fun, pars, fval);
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
    output.funcCount, fun, pars, fval);
end


