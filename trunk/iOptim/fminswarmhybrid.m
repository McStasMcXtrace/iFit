function [pars,fval,exitflag,output] = fminswarmhybrid(fun, pars, options,constraints, ub)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = FMINSWARMHYBRID(FUN,PARS,[OPTIONS],[CONSTRAINTS]) hybrid Particle Swarm Optimization
%
% This minimization method uses a hybrid Particle Swarm Optimization algorithm for 
% finding the minimum of the function 'FUN' in the real space. At each iteration 
% step, a local optimization is performed.
% Default local optimizer is the Nelder-Mead simplex (fminsearch). You may change
% it by defining the options.Hybrid function to any minimizer.
%
% Calling:
%   fminswarmhybrid(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminswarmhybrid(fun, pars, options) same as above, with customized options (optimset)
%   fminswarmhybrid(fun, pars, options, fixed) 
%     is used to fix some of the parameters. The 'fixed' vector is then 0 for
%     free parameters, and 1 otherwise.
%   fminswarmhybrid(fun, pars, options, lb, ub) 
%     is used to set the minimal and maximal parameter bounds, as vectors.
%   fminswarmhybrid(fun, pars, options, constraints) 
%     where constraints is a structure (see below).
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminswarmhybrid(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%      o=fminswarmhybrid('defaults');
%   options.Hybrid specifies the algorithm to use for local hybrid optimizations.
%      It may be set to any optimization method using the @fminsearch syntax.
%   option.SwarmSize sets the number of particules in the swarm (20-40).
%   option.SwarmC1 sets the local attractors strength (1-3)
%   option.SwarmC2 sets the global attractor strength (1-3).
%   option.SwarmW  sets inertia weight (0-1).
%
%  CONSTRAINTS may be specified as a structure
%   constraints.min= vector of minimal values for parameters
%   constraints.max= vector of maximal values for parameters
%   constraints.fixed= vector having 0 where parameters are free, 1 otherwise
%   constraints.step=  vector of maximal parameter changes per iteration
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% Reference:
% Kennedy J., Eberhart R.C. (1995): Particle swarm optimization. In: Proc.
% IEEE Conf. on Neural Networks, IV, Piscataway, NJ, pp. 1942-1948
% Shi, Y. and Eberhart, R. C. A modified particle swarm optimizer. Proc. 
% IEEE Int Conf. on Evol Comput pp. 69-73. IEEE Press, Piscataway, NJ, 1998
%
% Contrib:
% Alexandros Leontitsis leoaleq@yahoo.com Ioannina, Greece 2004
% and more informations on http://www.particleswarm.net, http://www.swarmintelligence.org
%
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='off';
  options.TolFun =1e-3;
  options.TolX   =1e-5;
  options.MaxIter=400;
  options.MaxFunEvals=400*50;
  options.Hybrid = @fminsearch;
  options.SwarmC1=2;
  options.SwarmC2=2;
  options.SwarmW =0;
  options.SwarmSize=20;
  pars = options;
  return
end

if nargin <= 2
	options=[];
end
if isempty(options)
  options=feval(mfilename, 'defaults');
end
if ~isfield(options,'Hybrid'), options.Hybrid=''; end
if isempty(options.Hybrid),    options.Hybrid=@fminsearch; end
if ~isfield(options,'SwarmC1'),options.SwarmC1=2; end
if ~isfield(options,'SwarmC2'),options.SwarmC2=2; end
if ~isfield(options,'SwarmW'), options.SwarmW =0; end
if ~isfield(options,'SwarmSize'), options.SwarmSize=20; end
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
end
if isfield(constraints, 'max')  % test if max values are valid
  index=find(isnan(constraints.max) | isinf(constraints.min));
  constraints.max(index) = 2*abs(pars(index));
end
if ~isfield(constraints, 'min')
  constraints.min = -2*abs(pars); % default min values
end
if ~isfield(constraints, 'max')
  constraints.max =  2*abs(pars); % default max values
end
if isfield(constraints, 'fixed') % fix some of the parameters if requested
  index = find(fixed);
  constraints.min(index) = pars(index); 
  constraints.max(index) = pars(index); % fix some parameters
end
options
% transfer optimset options and constraints
hoptions.space     = [ constraints.min(:) constraints.max(:) ];
hoptions.MaxIter   = options.MaxIter;
hoptions.Hybrid    = options.Hybrid;
hoptions.TolFun    = options.TolFun;
hoptions.Display   = options.Display;
hoptions.MaxFunEvals=options.MaxFunEvals;
hoptions.FunValCheck=options.FunValCheck;
hoptions.OutputFcn  =options.OutputFcn;
hoptions.c1         =options.SwarmC1;
hoptions.c2         =options.SwarmC2;
hoptions.w          =options.SwarmW;
hoptions.bees       =options.SwarmSize;
if isa(options.Hybrid, 'function_handle') | exist(options.Hybrid) == 2
  hoptions.StallFliLimit = 50;
else
  hoptions.StallFliLimit = Inf;
end
if isfield(constraints,'step')
  hoptions.maxv = constraints.step;
else
  hoptions.maxv = abs(constraints.max(:)-constraints.min(:))/2;
end

% call the optimizer
[pars,fval,exitflag,output] = hPSO(fun, pars, hoptions);


function [x,fval,istop,output]=hPSO(fitnessfun,pars,options,varargin)
%Syntax: [x,fval,exitflag,output]=hPSO(fitnessfun,nvars,options,varargin)
%___________________________________________________________________
%
% A hybrid Particle Swarm Optimization algorithm for finding the minimum of
% the function 'fitnessfun' in the real space.
%
% x is the scalar/vector of the functon minimum
% fval is the function minimum
% gfx contains the best particle for each flight (columns 2:end) and the
%  corresponding solutions (1st column)
% output structure contains the following information
%   reason : is the reason for stopping
%   flights: the nuber of flights before stopping
%   time   : the total time before stopping
% fitnessfun is the function to be minimized
% nvars is the number of variables of the problem
% options are specified in the file "PSOoptions.m"
%
%
% Reference:
% Kennedy J., Eberhart R.C. (1995): Particle swarm optimization. In: Proc.
% IEEE Conf. on Neural Networks, IV, Piscataway, NJ, pp. 1942-1948
%
%
% Alexandros Leontitsis
% Department of Education
% University of Ioannina
% 45110- Dourouti
% Ioannina
% Greece
% 
% University e-mail: me00743@cc.uoi.gr
% Lifetime e-mail: leoaleq@yahoo.com
% Homepage: http://www.geocities.com/CapeCanaveral/Lab/1421
% 
% 17 Nov, 2004.
nvars=length(pars);

if size(options.space,1)==1
    options.space=kron(options.space(:)', ones(nvars,1));
elseif size(options.space,1)~=nvars
    error('The rows of options.space are not equal to nvars.');
end

if size(options.maxv,1)==1
    options.maxv=maxv*ones(nvars,1);
elseif size(options.maxv,1)~=nvars
    error('The rows of options.maxv are not equal to nvars.');
end

c1    = options.c1;
c2    = options.c2;
w     = options.w;
maxv  = options.maxv;
space = options.space;
popul = options.bees;
flights = options.MaxIter;
Goal  = options.TolFun;

funcount=0; istop=0;

% Initial population (random start)
ru=rand(popul,size(space,1));
pop=ones(popul,1)*space(:,1)'+ru.*(ones(popul,1)*(space(:,2)-space(:,1))');
%pop(1,:)=pars(:); % force first particule to be the starting parameters

% Hill climb of each solution (bee)
if isa(options.Hybrid, 'function_handle') | exist(options.Hybrid) == 2
  for i=1:popul
      [pop(i,:),fxi(i,1),dummy,out]=feval(options.Hybrid,fitnessfun,pop(i,:));
      funcount = funcount+out.funcCount;
  end
end
% constrains within search domain
pop=min(pop,ones(popul,1)*space(:,2)');
pop=max(pop,ones(popul,1)*space(:,1)');
for i=1:popul
  fxi(i,1)=feval(fitnessfun,pop(i,:));
end
funcount = funcount+1;

if strcmp(options.Display,'iter')
  fmin_private_disp_start(mfilename, fitnessfun, pars);
end

% Local minima
p=pop;
fxip=fxi;

% Initialize the velocities
v=zeros(popul,size(space,1));

% Isolate the best solution
[Y,I]=min(fxi);
gfx(1,:)=[Y pop(I,:)];
P=ones(popul,1)*pop(I,:);

StallFli = 0;
message = 'Optimization terminated: maximum number of flights reached.';

% For each flight
for i=2:flights
    
    % Estimate the velocities
    r1=rand(popul,size(space,1));
    r2=rand(popul,size(space,1));
    v=v*w+c1*r1.*(p-pop)+c2*r2.*(P-pop);
    v=max(v,-ones(popul,1)*maxv');
    v=min(v,ones(popul,1)*maxv');
    
    % Add the velocities to the population 
    pop=pop+v;
    
    % Drag the particles into the search space
    pop=min(pop,ones(popul,1)*space(:,2)');
    pop=max(pop,ones(popul,1)*space(:,1)');
    
    % Hill climb search for the new population
    pnew=p;
    fxipnew=fxip;
    if isa(options.Hybrid, 'function_handle') | exist(options.Hybrid) == 2
      for j=1:popul
          [pop(j,:),fxi(j,1),dummy,out]     =feval(options.Hybrid,fitnessfun,pop(j,:));
          funcount = funcount+out.funcCount;
          [pnew(j,:),fxipnew(j,1),dummy,out]=feval(options.Hybrid,fitnessfun,p(j,:));
          funcount = funcount+out.funcCount;
      end
    else
      pnew=p;
    end
    pop=min(pop,ones(popul,1)*space(:,2)');
    pop=max(pop,ones(popul,1)*space(:,1)');
    pnew=min(pnew,ones(popul,1)*space(:,2)');
    pnew=max(pnew,ones(popul,1)*space(:,1)');
    for j=1:popul
      fxi(j,:)    =feval(fitnessfun,pop(j,:));
      fxipnew(j,:)=feval(fitnessfun,p(j,:));
    end
    funcount = funcount+2*popul;
    
    % Min(fxi,fxip)
    s=find(fxi<fxip);
    p(s,:)=pop(s,:);
    fxip(s)=fxi(s);
    
    % Min(fxipnew,fxip);
    s=find(fxipnew<fxip);
    p(s,:)=pnew(s,:);
    fxip(s)=fxipnew(s);
    
    % Isolate the best solution
    [Y,I]=min(fxip);
    gfx(i,:)=[Y p(I,:)];
    P=ones(popul,1)*p(I,:);
    
    % std stopping conditions
    % Get the point that correspond to the minimum of the function
    x=gfx(end,2:end);
    % Get the minimum of the function
    fval=gfx(end,1);
    [istop, message] = fmin_private_std_check(x, fval, i, funcount, options);
    if strcmp(options.Display, 'iter')
      fmin_private_disp_iter(i, funcount, fitnessfun, x, fval);
    end
  
    if istop
      break
    end
    
    % Termination conditions
    if gfx(i,1)==gfx(i-1,1)
        StallFli = StallFli+1;
    end    
    if StallFli >= options.StallFliLimit
        message = 'Optimization terminated: Stall Flights Limit reached.';
        istop=-8;
        break;
    end
end
if istop==0, message='Algorithm terminated normally'; end
output.iterations= i;
output.message   = message;
output.funcCount = funcount;
if isa(options.Hybrid, 'function_handle') | exist(options.Hybrid) == 2
output.algorithm = [ 'hybrid Particule Swarm Optimizer [' mfilename '/' localChar(options.Hybrid) ']' ];
else
output.algorithm = [ 'Particule Swarm Optimizer [fminswarm]' ];
end

if (istop & strcmp(options.Display,'notify')) | ...
   strcmp(options.Display,'final') | strcmp(options.Display,'iter')
  fmin_private_disp_final(output.algorithm, output.message, output.iterations, ...
    output.funcCount, fitnessfun, pars, fval);
end


