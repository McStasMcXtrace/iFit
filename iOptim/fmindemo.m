function abstract=fmindemo(dim, verbose, repeat_num)
% FMINDEMO optimization cross-comparison
% 
% A systematic test of all optimization methods is performed using a set of test
% problems. Starting configurations are chosen randomly.
% The verbose flag when set displays individual detailed optimization results.
% The test can be repeated iteratively so that a Monte Carlo sampling of starting
% parameters give a better stastistical estimate of each method efficiency.
%
% Calling:
%   fmindemo(dimensionality=vector, verbose=0|1, repetitions)
% Example:
%   fmindemo([ 2 10 50],1)  % detailed test whth 2,5 and 10 parameters
%   fmindemo(2,[], 10)      % 2 parameters optimization repeated 10 times
%
% Contrib:
% Test functions from Nikolaus Hansen, 2001-2007. e-mail: hansen@bionik.tu-berlin.de

if nargin == 0
  abstract{    1} = fmindemo(1, '',50); 
  abstract{end+1} = fmindemo(2, '',50);
  abstract{end+1} = fmindemo(4, '',50); 
  abstract{end+1} = fmindemo(8, '',40); 
  abstract{end+1} = fmindemo(12,'',35);
  abstract{end+1} = fmindemo(16,'',30); 
  abstract{end+1} = fmindemo(32,'',20); 
  abstract{end+1} = fmindemo(48,'',10); 
  abstract{end+1} = fmindemo(64,'',10); 
  return
end

if nargin == 1
  if isstruct(dim)
    abstract = dim;
    % display for given dimensionality an abstract of tests
    % matrix: rows=tests columns=optimizers, values=sorting index
    disp(['Test results: success probability and solving time. Dimensionality=' num2str(abstract.dim) ]);
    for j=1:length(abstract.optimizers)
      opt=abstract.optimizers{j};
      opt=opt(5:end);
      fprintf(1, '%10s ', opt(1:min(10, length(opt)))); 
      fprintf(1, '%6.2g ',  abstract.score(j).*100./abstract.repeat_num/length(abstract.problems)*2);
      fprintf(1, '%6.2g ',  abstract.duration(j)/abstract.score(j));
      fprintf(1, '\n');
    end
    
    return
  else
    verbose=[];
  end
end
if isempty(verbose), verbose=0; end
if nargin <= 2
  repeat_num=[];
end
if isempty(repeat_num), repeat_num=1; end

if isnumeric(dim) & length(dim) > 1
  for i=1:length(dim)
    fmindemo(dim(i), verbose);
  end
  return
end

optimizers={ ...
    'fminpowell','fminralg','fmincmaes','fminsimplex','fminnewton',...
    'fmingradrand', 'fminimfil', 'fminbfgs', 'fminlm', ...
    'fminhooke', 'fminanneal', 'fminsimpsa','fminsce','fminpso',...
    'fminga', ...
    'fminkalman',...
    'fminsearchOS', 'fminsearch', ...
    'fminswarm','fminswarmhybrid'};

if dim <=10
  optimizers{end+1} = 'fminrand';  % this one is too slow
  optimizers{end+1} = 'fminmulti'; % this one is not accurate
end

problems={...       % Unimodal functions
  'fjens1',       1, ...
  'fsphere',      1, ...
  'fssphere',     0.5, ...
  'fspherenoise', 1, ...
  'fmixranks',    1, ...
  'fsphereoneax', 1, ...
  'frandsphere',  1, ...
  'fspherelb0',   1, ...
  'fspherehull',  500, ...
  'fellilb0',     1, ...
  'fcornersphere',1, ...
  'fsectorsphere',1, ...
  'fstepsphere',  1, ...
  'fstep',        5.12, ...
  'flnorm',       1, ...
  'fneumaier3',   5, ...
  'fchangingsphere',1, ...
  'flogsphere',   1, ...
  'fexpsphere',   1, ...
  'fbaluja',      0.16, ...
  'fschwefel',    1, ...
  'fcigar',       1, ...
  'fcigtab',      1, ...
  'ftablet',      1, ...
  'felli',        1, ...
  'fellitest',    1, ...
  'fellii',       1, ...
  'fplane',       1, ...
  'ftwoaxes',     1, ...
  'fparabR',      1, ...
  'fsharpR',      1, ...
  'frosen',       1, ...
  'frosenlin',    1, ...
  'frosenmodif',  1, ...
  'fschwefelrosen1', 10, ...
  'fschwefelrosen2', 10, ...
  'fdiffpow',     1, ...
  'fabsprod',     1, ...
  'ffloor',       1, ...
  'fmax',         1, ...
  'fbirastrigin', 1, ...
  'fackley',      32000, ...    % Multimodal functions
  'fbohachevsky', 15, ...
  'fconcentric',  600, ...
  'fgriewank',    600, ...
  'fspallpseudorastrigin',30, ...
  'frastrigin',   1, ...
  'fschaffer',    100, ...
  'fschwefelmult',500, ...
  'ftwomax',      5, ...
  'ftwomaxtwo',   10, ...
  'frand',        1};

if dim == 1
  fig=figure;
  % plot 1D parameter space for each function
  n=floor(sqrt(length(problems)));
  m=ceil(length(problems)/n);
  for index_problem=1:2:length(problems)  % loop on functions
    psize=problems{index_problem+1};
    y=zeros(100);
    for i1=1:100
      p1=-psize+(i1/100)*2*psize;
      f = feval(str2func(problems{index_problem}), [ p1 ]);
      if length(f(:)) ~= 1
        error([ problems{index_problem} ' returns results of wrong dimension ' num2str(size(f)) ]);
      end
      y(i1) = f;
    end
    fprintf(1, '%s min=%g max=%g\n', problems{index_problem}, min(min(y)), max(max(y)));
    %subplot(m,n,index_problem);
    figure(fig);
    plot(y); title(problems{index_problem});
    print('-depsc', [ problems{index_problem} '_' num2str(dim) ]);
  end
  close(fig)
elseif dim == 2
  % plot 2D parameter space for each function
  fig=figure;
  n=floor(sqrt(length(problems)));
  m=ceil(length(problems)/n);
  for index_problem=1:2:length(problems)  % loop on functions
    psize=problems{index_problem+1};
    y=zeros(100,100);
    for i1=1:100
      p1=-psize+(i1/100)*2*psize;
      for i2=1:100
        p2=-psize+(i2/100)*2*psize;
        f = feval(str2func(problems{index_problem}), [ p1 p2 ]);
        if length(f(:)) ~= 1
          error([ problems{index_problem} ' returns results of wrong dimension ' num2str(size(f)) ]);
        end
        y(i1,i2) = f;
      end
    end
    fprintf(1, '%s min=%g max=%g\n', problems{index_problem}, min(min(y)), max(max(y)));
    %subplot(m,n,index_problem);
    figure(fig);
    surf(y); title(problems{index_problem});
    print('-depsc', [ problems{index_problem} '_' num2str(dim) ]);
  end
  close(fig)
end

dimensionality=randn(1,dim);

w=warning;
warning off
fprintf(1,'\n');
disp([ 'Starting optimization demo with ' num2str(dim) ' parameters to find.' ]);
disp(  'id name MaxFun MaxIter TolFun Description');
for index=1:length(optimizers)
  options=feval(optimizers{index},'defaults');
  opt=optimizers{index};
  if isfield(options, 'algorithm')
    alg=options.algorithm;
  else
    alg=opt;
  end
  options=feval(optimizers{index},'defaults');
  if length(opt) > 5  % shorten function mane for display
    opt=opt(5:end);
    opt=opt(1:min(5, length(opt)));
  end
  options.MaxIter=250*dim;
  options.MaxFunEvals=2500*dim;
  options.TolX=0;
  options.TolFun=1e-4;
  fprintf(1, '%2i %5s %6i %6i %6.2g %s\n', index, opt, options.MaxFunEvals, options.MaxIter, options.TolFun, alg);
end
fprintf(1,'\n');

t1 = clock;
abstract.score    = zeros(1,length(optimizers));
abstract.duration = zeros(1,length(optimizers));
abstract.success  = 0;
abstract.problems = problems;

for index_repeat=1:repeat_num

  if repeat_num>1, disp([ 'Repetition ' num2str(index_repeat) '/' num2str(repeat_num) ]); end

  for index_problem=1:2:length(problems)  % loop on functions
    fprintf(1, 'Problem: %15s [%i parameters] ', ...
      problems{index_problem}, length(dimensionality) );
    %    disp(startpars)
    startpars = dimensionality*problems{index_problem+1};
    fun=str2func(problems{index_problem});
    for index=1:length(optimizers)  % loop on optimizers
      if ~verbose, fprintf(1,'%i', index); end
      t0 = clock;
      options=feval(optimizers{index},'defaults');
      %if verbose, options.Display='iter'; end
      maxit = options.MaxIter; if ischar(maxit), maxit=eval(maxit,'0'); end
      if isinf(maxit), maxit=0; end
      options.MaxIter=min(2000, max(250*dim, maxit));
      maxfn = options.MaxFunEvals; if ischar(maxfn), maxfn=eval(maxfn,'0'); end
      if isinf(maxfn), maxfn=0; end
      options.MaxFunEvals=min(20000, max(2500*dim, maxfn));
      options.TolX=0;
      try
        [pars, fval, flag, out] = feval(optimizers{index}, fun, startpars(:), options);
      catch
        pars=[]; fval=Inf; out.funCount=Inf; flag=-4;
        out.iterations=Inf; out.algorithm=[ 'ERROR: ' optimizers{index} ];
        lasterr
      end
      duration(index)=etime(clock,t0);
      funcount(index)=out.funcCount;
      iterate(index) =out.iterations;
      if isnan(fval) | isinf(fval), fval=Inf; flag=-4; end
      if fval <= options.TolFun, flag=-1; f='CFUN';
      elseif ~isinf(fval) & iterate(index)>=options.MaxIter, f='nITR';
      elseif ~isinf(fval) & funcount(index)>=options.MaxFunEvals, f='nFUN';
      else
        switch flag
        case -1,f='Cfun';
        case -2,f='nitr';
        case -3,f='nFUN';
        case -4,f='*INF';
        case -5,f='Cx  ';
        case {-12, -9, -8}, f='c';
        otherwise
          if fval < options.TolFun, f='C';
          else f='n'; end
        end
      end
      fvals(index)   =fval;
      flags{index}   =f;
      algo{index}    =out.algorithm;
      params{index}=pars;
      if ~verbose, fprintf(1,'%s ', f(1)); end
    end % for optimizers
    disp(' ');
    % sort results
    criteria=duration.*fvals.*fvals;
    index=strmatch('C', flags);
    criteria(index) = criteria(index)*100; % put converged methods first
    index=strmatch('c', flags);
    criteria(index) = criteria(index)*100; % put converged methods first
    [dummy, sorti] = sort(criteria);
    
    % compute median duration
    mtime = mean(duration)+std(duration);
    score=1:length(optimizers);
    for index=1:length(optimizers)
      f = flags{index};
      if f(1) == 'C' || f(1)=='c',   
        if duration(index) > mtime, f='-'; 
        else f='+'; end
      elseif f(1)=='n', f=' '; 
      else f='*'; end
      % store results into array
      if f=='+' | f=='-', score(index)=1;
      else score(index)=0; end
    end
    abstract.score    = abstract.score   +score;
    abstract.duration = abstract.duration+duration.*score;

  end % end problems
  abstract.repeat_num = repeat_num;
  abstract.dim        = dim;
  abstract.optimizers = optimizers;
end % repeat

% display results
fmindemo(abstract);

warning(w);

return

% ------------------------------------------------------------------------------

% TEST functions taken from CMA-ES

function f=fjens1(x)
%
% use population size about 2*N
%
  f = sum((x>0) .* x.^1, 1);
  if any(any(x<0))
    idx = sum(x < 0, 1) > 0;
    f(idx) = 1e3;
%    f = f + 1e3 * sum(x<0, 1);
%    f = f + 10 * sum((x<0) .* x.^2, 1);
    f(idx) = f(idx) + 1e-3*abs(randn(1,sum(idx)));
%    f(idx) = NaN;
  end
  f=abs(sum(f));

function f=fsphere(x)
  x = x(:);
  f = sum(x.^2,1);
  f=abs(f);

function f=fssphere(x)
  x = x(:);
  f=sqrt(sum(x.^2, 1));

%  lb = -0.512; ub = 512; 
%  xfeas = x; 
%  xfeas(x<lb) = lb;
%  xfeas(x>ub) = ub; 
%  f=sum(xfeas.^2, 1);
%  f = f + 1e-9 * sum((xfeas-x).^2); 
  f=abs(f);
  
function f=fspherenoise(x)
  x = x(:);
  if length(x) < 2, x = [ x ; 0]; end
  [N,popsi] = size(x);
%  x = x .* (1 +  0.3e-0 * randn(N, popsi)/(2*N)); % actuator noise
  fsum = 10.^(0*(0:N-1)/(N-1)) * x.^2; 
  f = 0*rand(1,1) ...
      + fsum ...
      + fsum .* (2*randn(1,popsi) ./ randn(1,popsi).^0 / (2*N)) ...
      + 1*fsum^0.9 .* 2*randn(1,popsi) / (2*N); % 
%  f = fsum; 
  f=abs(f);

function f=fmixranks(x)
  x=x(:);
  if length(x) < 2, x = [ x ; 0]; end
  N = length(x);
  f=(10.^(0*(0:(N-1))/(N-1))*x.^2).^0.5;
  if size(x, 2) > 1 % compute ranks, if it is a population 
    [ignore, idx] = sort(f);
    [ignore, ranks] = sort(idx);
    k = 9; % number of solutions randomly permuted, lambda/2-1
           % works still quite well (two time slower)
    for i = k+1:k-0:size(x,2)
      idx(i-k+(1:k)) = idx(i-k+randperm(k)); 
    end
    %disp([ranks' f'])
    [ignore, ranks] = sort(idx);
    %disp([ranks' f'])
    %pause
    f = ranks+1e-9*randn(1,1);
  end
  f=abs(f);
  
function f = fsphereoneax(x)
  x=x(:);
  f = x(1)^2;
  f = mean(x)^2;
  f=abs(f);
  
function f=frandsphere(x)
  x=x(:);
  N = length(x);
  idx = ceil(N*rand(7,1));
  f=sum(x(idx).^2);
  f=abs(f);

function f=fspherelb0(x, M) % lbound at zero for 1:M needed
  x=x(:);
  if nargin < 2 M = 0; end
  N = length(x);
  % M active bounds, f_i = 1 for x = 0
  f = -M + sum((x(1:M) + 1).^2);
  f = f + sum(x(M+1:N).^2);
  f=abs(f);
  
function f=fspherehull(x)
  x=x(:);
  % Patton, Dexter, Goodman, Punch
  % in -500..500
  % spherical ridge through zeros(N,1)
  % worst case start point seems x = 2*100*sqrt(N)
  % and small step size
  N = length(x);
  f = norm(x) + (norm(x-100*sqrt(N)) - 100*N)^2;
  f=abs(f);
  
function f=fellilb0(x, idxM, scal) % lbound at zero for 1:M needed
  x=x(:);
  if length(x) < 2, x = [ x ; 0]; end
  N = length(x);
  if nargin < 3 || isempty(scal)
    scal = 100;
  end
  scale=scal.^((0:N-1)/(N-1));
  if nargin < 2 || isempty(idxM)
    idxM = 1:N;
  end
  %scale(N) = 1e0;
  % M active bounds
  xopt = 0.1;
  x(idxM) = x(idxM) + xopt;
  f = scale.^2*x.^2;
  f = f - sum((xopt*scale(idxM)).^2); 
%  f = exp(f) - 1;
%  f = log10(f+1e-19) + 19;

  f = f + 1e-19;
  f=abs(f);
  
function f=fcornersphere(x)
  x=x(:);
  w = ones(length(x));
  w(1) = 2.5; w(2)=2.5;
  idx = x < 0;
  f = sum(x(idx).^2);
  idx = x > 0;
  f = f + 2^2*sum(w(idx).*x(idx).^2);
  f=abs(f);
  
function f=fsectorsphere(x, scal)
%
% This is deceptive for cumulative sigma control CSA in large dimension:
% The strategy (initially) diverges for N=50 and popsize = 150.  (Even
% for cs==1 this can be observed for larger settings of N and
% popsize.) The reason is obvious from the function topology. 
% Divergence can be avoided by setting boundaries or adding a
% penalty for large ||x||. Then, convergence can be observed again. 
% Conclusion: for popsize>N cumulative sigma control is not completely
% reasonable, but I do not know better alternatives. In particular:
% TPA takes longer to converge than CSA when the latter still works. 
%
  x=x(:);
  if nargin < 2 || isempty (scal)
    scal = 1e3;
  end
  f=sum(x.^2,1);
  idx = find(x<0);
  if ~isempty(idx)
    f = f + (scal-1)^2 * sum(x(idx).^2,1);
  end
  f=abs(f);
  
function f=fstepsphere(x, scal)
  x=x(:);
  if length(x) < 2, x = [ x ; 0]; end
  if nargin < 2 || isempty (scal)
    scal = 1e0;
  end
  N = length(x);
  f=1e-11+sum(scal.^((0:N-1)/(N-1))*floor(x+0.5).^2);
  f=1e-11+sum(floor(scal.^((0:N-1)/(N-1))'.*x+0.5).^2);
%  f=1e-11+sum(floor(x+0.5).^2);
  f=abs(f);

function f=fstep(x)
  x=x(:);
  % in -5.12..5.12 (bounded)
  N = length(x);
  f=1e-11+6*N+sum(floor(x));
  f=abs(f);

function f=flnorm(x, scal, e)
  x=x(:);
  if length(x) < 2, x = [ x ; 0]; end
  if nargin < 2 || isempty(scal)
    scal = 1;
  end
  if nargin < 3 || isempty(e)
    e = 1;
  end
  if e==inf
    f = max(abs(x));
  else
    N = length(x);
    scale = scal.^((0:N-1)/(N-1))';
    f=sum(abs(scale.*x).^e);
  end
  f=abs(f);

function f=fneumaier3(x) 
  x=x(:);
  % in -n^2..n^2
  % x^*-i = i(n+1-i)
  N = length(x);
%  f = N*(N+4)*(N-1)/6 + sum((x-1).^2) - sum(x(1:N-1).*x(2:N));
  f = sum((x-1).^2) - sum(x(1:N-1).*x(2:N));
  f=abs(f);
  
function f=fchangingsphere(x)
  x=x(:);
  N = length(x);
  global scale_G; global count_G; if isempty(count_G) count_G=-1; end
  count_G = count_G+1;
  if mod(count_G,10) == 0
    scale_G = 10.^(2*rand(1,N));
  end
  %disp(scale(1));
  f = scale_G*x.^2;
  f=abs(f);
  
function f= flogsphere(x)
  x=x(:);
  f = 1-exp(-sum(x.^2));
  
function f= fexpsphere(x)
  f = exp(sum(x.^2)) - 1;
  f=abs(f);
  
function f=fbaluja(x)
  % in [-0.16 0.16]
  x=x(:);
  y = x(1);
  for i = 2:length(x)
    y(i) = x(i) + y(i-1);
  end
  f = 1e5 - 1/(1e-5 + sum(abs(y)));
  f=abs(f);

function f=fschwefel(x)
  x=x(:);
  f = 0;
  for i = 1:length(x),
    f = f+sum(x(1:i))^2;
  end
  f=abs(f);

function f=fcigar(x, ar)
  x=x(:);
  if nargin < 2 || isempty(ar)
    ar = 1e3;
  end
  f = x(1,:).^2 + ar^2*sum(x(2:end,:).^2,1);
  f=abs(f);
  
function f=fcigtab(x)
  x=x(:);
  f = x(1,:).^2 + 1e8*x(end,:).^2 + 1e4*sum(x(2:(end-1),:).^2, 1);
  f=abs(f);
  
function f=ftablet(x)
  x=x(:);
  f = 1e6*x(1,:).^2 + sum(x(2:end,:).^2, 1);
  f=abs(f);

function f=felli(x, lgscal, expon, expon2)
  x=x(:);
  % lgscal: log10(axis ratio)
  % expon: x_i^expon, sphere==2
  if length(x) < 2, x = [ x ; 0]; end
  N = length(x);

%  x = x - repmat(-0.5+(1:N)',1,size(x,2)); % optimum in 1:N
  if nargin < 2 || isempty(lgscal), lgscal = 3; end
  if nargin < 3 || isempty(expon), expon = 2; end
  if nargin < 4 || isempty(expon2), expon2 = 1; end

  f=((10^(lgscal*expon)).^((0:N-1)/(N-1)) * abs(x).^expon).^(1/expon2);
%  if rand(1,1) > 0.015
%    f = NaN;
%  end
%  f = f + randn(size(f));
  f=abs(f);

function f=fellitest(x)
  x=x(:);
  if length(x) < 2, x = [ x ; 0]; end
  beta = 0.9;
  N = length(x);
  f = (1e6.^((0:(N-1))/(N-1))).^beta * (x.^2).^beta; 
  f=abs(f);
  
function f=fellii(x, scal)
  x=x(:);
  if length(x) < 2, x = [ x ; 0]; end
  N = length(x);
  if nargin < 2
    scal = 1;
  end
  f= (scal*(1:N)).^2 * (x).^2;
  f=abs(f);

function f=fplane(x)
  x=x(:);
  f=x(1);
  f=abs(f);

function f=ftwoaxes(x)
  x=x(:);
  f = sum(x(1:floor(end/2)).^2) + 1e6*sum(x(floor(1+end/2):end).^2);
  f=abs(f);

function f=fparabR(x)
  x=x(:);
  f = -x(1,:) + 100*sum(x(2:end,:).^2,1);
  f=abs(f);

function f=fsharpR(x)
  x=x(:);
  f = abs(-x(1)) + 100*norm(x(1:end));
  f=abs(f);
  
function f=frosen(x)
  x=x(:);
  if length(x) < 2, x = [ x ; 0]; end
  N = length(x); 
  popsi = size(x,2); 
  f = 1e2*sum((x(1:end-1,:).^2 - x(2:end,:)).^2,1) + sum((x(1:end-1,:)-1).^2,1);
  % f = f + f^0.9 .* (2*randn(1,popsi) ./ randn(1,popsi).^0 / (2*N)); 
  f=abs(f);

function f=frosenlin(x)
  x=x(:);
  if length(x) < 2, x = [ x ; 0]; end

  x_org = x;
  x(x>30) = 30;
  x(x<-30) = -30;

  f = 1e2*sum(-(x(1:end-1,:).^2 - x(2:end,:)),1) + ...
      sum((x(1:end-1,:)-1).^2,1);

  f = f + sum((x-x_org).^2,1);
%  f(any(abs(x)>30,1)) = NaN; 
  f=abs(f);

function f=frosenmodif(x)
  x=x(:);
  if length(x) < 2, x = [ x ; 0]; end
  f = 74 + 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 ...
      - 400*exp(-sum((x+1).^2)/2/0.05);
  f=abs(f);
  
function f=fschwefelrosen1(x)
  % in [-10 10] 
  x=x(:);
  f=sum((x.^2-x(1)).^2 + (x-1).^2);
  f=abs(f);
  
function f=fschwefelrosen2(x)
  % in [-10 10] 
  x=x(:);
  if length(x) < 2, x = [ x ; 0]; end
  f=sum((x(2:end).^2-x(1)).^2 + (x(2:end)-1).^2);
  f=abs(f);

function f=fdiffpow(x)
  x=x(:);
  if length(x) < 2, x = [ x ; 0]; end
  N = length(x); 
  f=sum(abs(x).^(2+10*(0:N-1)'/(N-1)));
  % f = sqrt(f); 
  f=abs(f);
  
function f=fabsprod(x)
  x=x(:);
  f = sum(abs(x),1) + prod(abs(x),1);
  f=abs(f);

function f=ffloor(x)
  x=x(:);
  f = sum(floor(x+0.5).^2,1); 
  f=abs(f);

function f=fmax(x)
  x=x(:);
  f = max(abs(x), [], 1);
  f=abs(f);

%%% Multimodal functions 

function f=fbirastrigin(x)
% todo: the volume needs to be a constant 
  x=x(:);
  N = length(x); 
  idx = (sum(x, 1) < 0.5*N); % global optimum
  f = zeros(1,size(x,2));
  f(idx) = 10*(N-sum(cos(2*pi*x(:,idx)),1)) + sum(x(:,idx).^2,1); 
  idx = ~idx;
  f(idx) = 0.1 + 10*(N-sum(cos(2*pi*(x(:,idx)-2)),1)) + sum((x(:,idx)-2).^2,1); 
  f=abs(f);

function f=fackley(x)
  x=x(:);
  % -32.768..32.768
  % Adding a penalty outside the interval is recommended,  
  % because for large step sizes, fackley imposes like frand
  % 
  N = length(x); 
  f = 20-20*exp(-0.2*sqrt(sum(x.^2)/N)); 
  f = f + (exp(1) - exp(sum(cos(2*pi*x))/N));
  % add penalty outside the search interval
  f = f + sum((x(x>32.768)-32.768).^2) + sum((x(x<-32.768)+32.768).^2);
  f=abs(f);
  
function f = fbohachevsky(x)
 % -15..15
  x=x(:);
  if length(x) < 2, x = [ x ; 0]; end
  f = sum(x(1:end-1).^2 + 2 * x(2:end).^2 - 0.3 * cos(3*pi*x(1:end-1)) ...
	  - 0.4 * cos(4*pi*x(2:end)) + 0.7);
	f=abs(f);
  
function f=fconcentric(x)
  % in  +-600
    x=x(:);
  s = sum(x.^2);
  f = s^0.25 * (sin(50*s^0.1)^2 + 1);
  f=abs(f);

function f=fgriewank(x)
  % in [-600 600]
  x=x(:);
  N = length(x);
  f = 1 - prod(cos(x'./sqrt(1:N))) + sum(x.^2)/4e3;
  % f = f + 1e4*sum(x(abs(x)>5).^2);
  % if sum(x(abs(x)>5).^2) > 0
  %   f = 1e4 * sum(x(abs(x)>5).^2) + 1e8 * sum(x(x>5)).^2;
  % end
  f=abs(f);
  
function f=fspallpseudorastrigin(x, scal, skewfac, skewstart, amplitude)
% by default multi-modal about between -30 and 30
  x=x(:);
  if nargin < 5 || isempty(amplitude)
    amplitude = 40;
  end
  if nargin < 4 || isempty(skewstart)
    skewstart = 0;
  end
  if nargin < 3 || isempty(skewfac)
    skewfac = 1;
  end
  if nargin < 2 || isempty(scal)
    scal = 1;
  end
  N = length(x); 
  scale = 1;
  if N > 1
    scale=scal.^((0:N-1)'/(N-1)); 
  end
  % simple version: 
  % f = amplitude*(N - sum(cos(2*pi*(scale.*x)))) + sum((scale.*x).^2);

  % skew version: 
  y = repmat(scale, 1, size(x,2)) .* x;
  idx = find(x > skewstart);
  if ~isempty(idx)
    y(idx) =  skewfac*y(idx);
  end
  f = amplitude * (0*N-prod(cos((2*pi)^0*y),1)) + 0.05 * sum(y.^2,1) ...
      + randn(1,1);
  f=abs(f);

function f=frastrigin(x, scal, skewfac, skewstart, amplitude)
% by default multi-modal about between -30 and 30
  x=x(:);
  if nargin < 5 || isempty(amplitude)
    amplitude = 10;
  end
  if nargin < 4 || isempty(skewstart)
    skewstart = 0;
  end
  if nargin < 3 || isempty(skewfac)
    skewfac = 1;
  end
  if nargin < 2 || isempty(scal)
    scal = 1;
  end
  N = length(x); 
  scale = 1;
  if N > 1
    scale=scal.^((0:N-1)'/(N-1)); 
  end
  % simple version: 
  % f = amplitude*(N - sum(cos(2*pi*(scale.*x)))) + sum((scale.*x).^2);

  % skew version: 
  y = repmat(scale, 1, size(x,2)) .* x;
  idx = find(x > skewstart);
  % idx = intersect(idx, 2:2:10); 
  if ~isempty(idx)
    y(idx) =  skewfac*y(idx);
  end
  f = amplitude * (N-sum(cos(2*pi*y),1)) + sum(y.^2,1);
  f=abs(f);
  
function f = fschaffer(x)
 % -100..100
  x=x(:);
  if length(x) < 2, x = [ x ; 0]; end
  N = length(x);
  s = x(1:N-1).^2 + x(2:N).^2;
  f = sum(s.^0.25 .* (sin(50*s.^0.1).^2+1));
  f=abs(f);

function f=fschwefelmult(x)
  % -500..500
  % 
  x=x(:);
  N = length(x); 
  f = - sum(x.*sin(sqrt(abs(x))), 1);
  % f = 418.9829*N - 1.27275661e-5*N - sum(x.*sin(sqrt(abs(x))), 1);
  % penalty term 
  y = 1e4*sum((abs(x(abs(x)>500))-500).^2, 1);
  if ~isempty(y)
    f = f + y;
  end
  f=abs(f);
  
function f=ftwomax(x)
  % Boundaries at +/-5
    x=x(:);
  N = length(x); 
  f = -abs(sum(x)) + 5*N;
  f=abs(f);

function f=ftwomaxtwo(x)
  % Boundaries at +/-10
    x=x(:);
  N = length(x); 
  f = abs(sum(x));
  if f > 30
    f = f - 30;
  end
  f = -f;
  f=abs(f);

function f=frand(x)
  f=1/(1-rand) - 1;
  f=abs(f);

  
function n=norm2(x)
x = x(:);
n=sqrt(sum(abs(x).*abs(x)));


