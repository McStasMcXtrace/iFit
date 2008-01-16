function [pars,fval,message,output] = fminga(fun, pars, options)
% [MINIMUM,FVAL,MESSAGE,OUTPUT] = FMINGA(FUN,PARS,[OPTIONS]) genetic algorithm optimizer
%
% This optimization method uses a Genetic Algorithm (real coding)
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
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          MESSAGE return state of the optimizer
%          OUTPUT additional information returned as a structure.
% Contrib:
% By: Javad Ivakpour
% E-mail: javad7@gmail.com
% May 2006
%
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='off';
  options.TolFun =1e-3;
  options.TolX   =1e-5;
  options.MaxIter=3000;
  options.MaxFunEvals=3000*100;
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
[pars,fval,message,output] = GA(fun, pars, options);

% private function ------------------------------------------------------------


function [pars, fval, message, output] = GA(fun, pars, options)
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
                  %   number of variables that is used in the function in
                  %   fun00.m file)
n=100;            % Number of population

m0=options.MaxIter;            % Number of generations that max value remains constant
                  %   (use for termination criteria)
nmutationG=20;                  %number of mutation children(Gaussian)
nmutationR=20;                  %number of mutation children(random)
nelit=2;                        %number of elitism children
%-------------------------------------------------------------------------
nmutation= nmutationG+nmutationR;
max1     = zeros(nelit,var);
parent   = zeros(n,var);
p        = zeros(n,var);
stdp     = 0.05;
sigma=pars.*stdp/5;    %Parameter that related to Gaussian
                        %   function and used in mutation step
for l=1:var
    p(:,l)=pars(l)*(1+stdp*randn(n,1));
end
p(1,:)=pars;
m=m0;
maxvalue=ones(m,1)*-1e10;
maxvalue(m)=-1e5;
g=0;
iterations=0;
funcount  =0;
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
    maxvalue1(1:nelit)=s(n:-1:n-nelit+1);
    if nelit==0
        maxvalue1(1)=s(n);
        for i=1:n
            if y(i)==maxvalue1(1)
                max1(1,:)=p(i,:);
            end
        end
    end
    for k=1:nelit
        for i=1:n
            if y(i)==maxvalue1(k)
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
            p(i,1:var)=p(i,1:var)+rand(1,var).*sigma;
        end
    end
    m=m+1;
    maxvalue(m)=maxvalue1(1);
    
    % std stopping conditions
    fval = -maxvalue1(1);
    pars =  max1(1,:); pars = pars(:)';
    iterations = iterations+1;
    funcount = n*iterations;
    [istop, message] = fmin_private_std_check(pars, fval, iterations, funcount, options);
  
    if istop
      break
    end
end

% output results --------------------------------------------------------------
output.iterations = iterations;
output.algorithm  = mfilename;
output.message    = message;
output.funcCount  = funcount;

