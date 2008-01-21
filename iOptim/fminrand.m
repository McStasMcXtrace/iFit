function [pars,fval,exitflag,output] = fminrand(fun, pars, options)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = FMINRAND(FUN,PARS,[OPTIONS]) adaptive random search optimizer
%
% This minimization method uses an adaptive random search.
% 
% Calling:
%   fminrand(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminrand(fun, pars, options) same as above, with customized options (optimset)
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
%   optimset('fminrand')
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% Contrib: Argimiro R. Secchi (arge@enq.ufrgs.br) 2001
% Modified by Giovani Tonel(giovani.tonel@ufrgs.br) on September 2006
%   A.R. Secchi and C.A. Perlingeiro, "Busca Aleatoria Adaptativa",
%   in Proc. of XII Congresso Nacional de Matematica Aplicada e
%   Computacional, Sao Jose do Rio Preto, SP, pp. 49-52 (1989).
%
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='off';
  options.TolFun =1e-6;
  options.TolX   =1e-5;
  options.MaxIter=500;
  options.MaxFunEvals=5000;
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
[pars,fval,exitflag,output] = buscarnd(fun, pars, options);

% private function ------------------------------------------------------------

function [xo,Ot,istop,output]=buscarnd(S,x0,options)
%   Unconstrained global optimization using adaptive random search.
%
%   [xo,Ot,nS]=buscarnd(S,x0,options.Display,nOt,samples,Lb,Ub,problem,options.TolX,options.MaxIter,R,red,mem)
%

%   Copyright (c) 2001 by LASIM-DEQUI-UFRGS
%   $Revision: 1.2 $  $Date: 2008-01-21 09:04:07 $
%   Argimiro R. Secchi (arge@enq.ufrgs.br)
%
%   Based on the algorithm of the same author written in C
%   for constrained global optimization in 1989/09/29 and published as:
%   A.R. Secchi and C.A. Perlingeiro, "Busca Aleatoria Adaptativa",
%   in Proc. of XII Congresso Nacional de Matematica Aplicada e
%   Computacional, Sao Jose do Rio Preto, SP, pp. 49-52 (1989).
%
% Modified by Giovani Tonel(giovani.tonel@ufrgs.br) on September 2006



rand('state',sum(100*clock));   %rand('state',sum(100*clock)) resets it to
                          %             a different state each time.

%   S: objective function
%   x0: initial point
%   options.Display: (0): no plot (default), (>0) plot figure options.Display with pause, (<0) plot figure options.Display
%   nOt: maximum number of optimal points (default = 1)
%   samples: number of samples per stage (default = 3*size(x0(:),1))
%   Lb, Ub: lower and upper bound vectors to plot (default = x0*(1+/-2))
%   problem: (-1): minimum (default), (1): maximum
%   options.TolX: tolerance (default = 1e-4)
%   options.MaxIter: maximum number of stages (default = 50*(1+4*~(options.Display>0)))
%   R: axis vector of the hyperellipse centered in x0 (default = max(0.1*abs(x0+~x0),1))
%   red: factor to reduce the size of the axis vector (0,1) (default = 0.2)
%   mem: number of stored points during the exploration phase [0, min(18,samples-1)]
%   xo: matrix of optimal points (column vectors)
%   Ot: vector of optimal values of S
%   nS: number of objective function evaluations

nOt=1;
samples=3*size(x0(:),1);
problem=-1; istop=0;
R=max(0.1*abs(x0+~x0),1);

red=0.2;

mem=min(18,samples-1);

mem=samples-1;

                    % initialization 
x0=x0(:);
R=abs(R(:));
n=size(x0,1);
distrib=1;          % distr =  1  --> exploration phase
                   % distr >  1  --> progression phase
asym1=-1*ones(n,1); % asym1 =  1  --> asymmetric distribution to reduce xo
                   % asym1 = -1  --> asymmetric distribution to increase xo
asym2=2*ones(n,1);  % asym2 =  2  --> symmetric distribution
                   % asym2 = 1.5 --> asymmetric distribution
metric=0.5*norm2(R);
a0=0.5*metric;
a2=0.1*metric;
a3=10*options.TolX;
n3=1;
holdm=floor(1/red+0.5)*n;
nhold=0;            % counter to keep distribution constant during holdm times
top=0;
if mem > 0,
  n4=10*mem;
  idx=[1:n4+1];
  last=n4;
  first=0;
  idx(last+1)=first;
  mem_x=zeros(n,n4);
  mem_S=-ones(1,n4)*inf;
end

%     Criterion of axis vector reduction
%
%     R_red = R*0.5^(distrib-1)
%
%     or norm2(R_red) = norm2(R)*0.5^(distrib-1)
%
%     then if the criterion is norm2(R_red) < options.TolX = 10^(-asc)
%
%     the distrib is limited by
%
%     distrib = 1 + (asc + log10(norm2(R))) / log10(2)

lim_distrib=floor(0.5*(1+log10(metric/options.TolX)/log10(2.0)));

xOt=[];
Ot=[];
xo=x0;
yo=feval(S,x0)*problem;
if strcmp(options.Display,'iter')
  fmin_private_disp_start(mfilename, S, xo, yo*problem);
end

nS=1;
it=0;
opt=0;
x=zeros(n,samples);
y=zeros(1,samples);
 
while 1 % opt < nOt %  it < options.MaxIter & opt < nOt,
  l3=0;
  it=it+1;
  for j=1:samples,   % sampling
    a=min(max(rand(n,1),0.05),0.95);
    x(:,j)=xo+(R.*asym1.*(asym2.*a-1).^distrib)/distrib;
    y(j)=feval(S,x(:,j))*problem;
    nS=nS+1;
  end

  l4=mem & distrib == 1;  % sorting the samples
  if samples > 1,
    if ~l4, 
      [y(1),i]=max(y);
      x(:,1)=x(:,i);
    else
      [y,i]=sort(-y); y=-y;
      x=x(:,i);
    end
  end

  if l4,                  % memorize samples
    n2=top;
    l1=0;
    for j=1:samples,
      n1=0;
      if opt,
        for i=1:opt
          if norm2(xOt(:,i)-x(:,j)) < a2,
            n1=1;
            if j == 1,
              l1=1;
            end
            break;
          end
        end
      end
      
      if ~n1,       
        for i=1:j-1
          if norm2(x(:,j)-x(:,i)) < metric,
            n1=1;
            break;
          end
        end
  
        if ~n1,
          k=idx(1);
          for i=1:top
            if norm2(mem_x(:,k)-x(:,j)) < metric,
              n1=1;
              break;
            else
              k=idx(k+1);
            end
          end
          
          if ~n1,
            first=idx(first+1);
            if ~first,              % expand memory
              mem_x=[mem_x,zeros(n,n4)];
              mem_S=[mem_S,-ones(1,n4)*inf];
              idx=[idx,[n2+1:n2+n4+1]];
              first=n2;
              idx(last+1)=first;
              last=n2+n4;
              idx(last+1)=0;
            end
            k=first;
            n2=n2+1;
          end  

          if ~n1 | y(j) > mem_S(k),
            mem_x(:,k)=x(:,j);
            mem_S(:,k)=y(j);
          end
          
          if l1,
            l1=0;
            x(:,1)=x(:,j);
            y(1)=y(j);
          end    
        end  % if ~n1
      end % if ~n1
          
      if n2-top == mem,
        break;
      end
    end % for j
    top=n2;
  end % if l4
                  % analysis of best sampled point
  a1=norm2(x(:,1)-xo)/(0.1+norm2(xo))+abs(y(1)-yo)/(0.1+abs(yo));
  l1=y(1) > yo;
  if l1,
    if norm2(x(:,1)-xo) > a0,
      nhold=0;
      n3=1;
    end
    
    z=xo; xo=x(:,1); x(:,1)=z;
    a4=yo; yo=y(1); y(1)=a4;
  end

  if a1 < options.TolX & distrib >= lim_distrib,
    l2=1;
    if opt,
      for i=1:opt,
        if norm2(xOt(:,i)-xo) < a3,
          l2=0;
          break;
        end
      end
    end

    if l2,
      opt=opt+1;
      Ot=[Ot,yo];
      xOt=[xOt,xo];
    end

        % rescue another point from memory
    if ~top,
      if opt,
        xo=xOt;
        yo=Ot;
      end
      % break;
    end
    
    if l2 & ~l1,
      n2=0;
      for j=1:top,
        k=idx(n2+1);
        if norm2(mem_x(:,k)-xo) < a2,
          top=top-1;
          idx(last+1)=k;
          last=k;
          idx(n2+1)=idx(k+1);
          idx(k+1)=0;
          if k == first,
            first=n2;
          end
        else
          n2=k;
        end
      end
    end

    if ~top,
      if opt,
        xo=xOt;
        yo=Ot;
      end
      % break;
    end  
          
    k=idx(1);
    xo=mem_x(:,k);
    yo=mem_S(k);
    idx(last+1)=k;
    last=k;
    idx(1)=idx(k+1);
    idx(k+1)=1;
    top=top-1;
    
    if k == first,
      first=1;
    end
    
    l3=1;
  elseif opt & (~l1 | (l1 & ~l4))
    n1=0;
    for i=1:opt,
      if norm2(xOt(:,i)-xo) < a2,
        n1=1;
        break;
      end
    end
  
    if n1,
      if ~top,
        xo=xOt;
        yo=Ot;
        % break;
      end

      k=idx(1);
      xo=mem_x(:,k);
      yo=mem_S(k);
      idx(last+1)=k;
      last=k;
      idx(1)=idx(k+1);
      idx(k+1)=1;
      top=top-1;

      if k == first,
        first=1;
      end

      l3=1;
    end
  end % if a1
                 % adjust search direction and criterion
  if l3,
    n3=1;
    distrib=1;
    nhold=0;
  elseif (a1 < a2 | ~l1) & (nhold+n3 > holdm),
    n3 = n3+~mod(distrib,3);
    distrib=distrib+2;
    nhold=0;
  else
    nhold=nhold+n3;
  end
  
  for i=1:n,
    if l3 | (~l3 & x(i,1) == xo(i)),
      asym2(i)=2;
    else
      asym2(i)=1.5;
      if x(i,1) < xo(i),
        asym1(i)=1;
      else
        asym1(i)=-1;
      end
    end
  end
  
  % std stopping conditions
  [istop, message] = fmin_private_std_check(xo, yo*problem, it, nS, options);
  if strcmp(options.Display, 'iter')
    fmin_private_disp_iter(it, nS, S, x, yo*problem);
  end
  
  if istop
    opt = nOt;
    break
  end
end % while
 
if opt == nOt,
  yo=Ot;
  xo=xOt;
end

% get the best solution
if opt >= 1,
  [yo,i]=sort(-yo);
  yo=-yo(1);
  xo=xo(:,i(1));
end

Ot=yo*problem;
xo = xo(:)';

pars=xo;
fval=Ot;

% output results --------------------------------------------------------------
if istop==0, message='Algorithm terminated normally'; end
output.iterations = it;
output.algorithm  = [ 'Adaptive Random Search [' mfilename ']' ];
output.message    = message;
output.funcCount  = nS;

if (istop & strcmp(options.Display,'notify')) | ...
   strcmp(options.Display,'final') | strcmp(options.Display,'iter')
  fmin_private_disp_final(output.algorithm, output.message, output.iterations, ...
    output.funcCount, fun, pars, fval);
end

function n=norm2(x)
x = x(:);
n=sqrt(sum(abs(x).*abs(x)));

